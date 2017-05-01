#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Take any two lists of SNPs as rsID or chr:location and calculate pairwise LD.

Pairwise LD is calculated between each SNP in the first list (SNP1) and every
SNP is the second list within n kilobases of the first SNP.

The list is filtered by an R\u00b2 greater than x and the set of all possible
SNPs is returned as a set of SNP_Link objects to allow lookup.

Alternatively, the links can be output to a pandas/written table in the format:
    SNP1, SNP2, R2, SNP1_A1, SNP1_A2, SNP2_A1, SNP2_A2, SNP1_A1_LINKED_WITH

The final column in the above will be either 'A1' or 'A2' and defines which
allele in SNP2 is linked with which allele in SNP1.

If run from the command line, either the tab-delimited or pandas format table
may be written to a file, no SNP_Link objects will be returned.

Note: these calculations only work on bi-allelic SNPs, and more complex SNPs
or indels are ignored.

Note: this code currently converts all rsIDs to chr:pos format and keeps them
that way, so rsIDs are lost during parsing and only positions are returned, to
get rsIDs back use dbSNP.lookup_locations().

WARNING: Given an average linkage of 10 SNPs in list 2 to every SNP in list 1,
         if list 1 is large, the output can be a huge memory hog. In this
         situation writing to a table is a better choice.
         Also, if SNPs are given as rsIDs, every SNP will be converted to a
         location by a database lookup. This is quite fast, but if both lists
         are large, it makes sense to only do the conversion once. Use the
         convert_rsid_lists() function for that.

Runtime
-------
The lookup is parallelized by (fyrd)[https://fyrd.science] if available, or
using local cores only.

If using fyrd, the initial SNP list conversion and filtration is done in
parallel on the local machine only, the lookups are then bundled into sets of
10000 jobs and submitted to the cluster. Given a list1 of 40,000 SNPs, 400,000
total lookups will be done, which will result in 4,000 jobs submitted to the
cluster, each will take about 90 minutes to run on average hardware

Examples
--------
TODO
"""
import os as _os
import re as _re
import sys as _sys
import argparse as _argparse
from tempfile import mkstemp as _temp

import dbSNP as _dbSNP
import subprocess as _sub

DB_PATH = '/godot/dbsnp'
DB_VERS = 150

POPULATIONS = ["ALL", "AFR", "AMR", "EAS", "EUR", "SAS", "ACB", "ASW", "BEB",
               "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH",
               "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL",
               "PJL", "PUR", "STU", "TSI", "YRI"]

# Set data directories
DATA_DIR = "/godot/1000genomes/1000GP_Phase3"


class MissingSNPError(Exception):

    """Exceception to catch missing SNPs."""

    pass


class BadSNPError(Exception):

    """Exceception to catch missing SNPs."""

    pass


def run(command, raise_on_error=False):
    """Run a command with subprocess the way it should be.

    Parameters
    ----------
    command : str
        A command to execute, piping is fine.
    raise_on_error : bool
        Raise a subprocess.CalledProcessError on exit_code != 0

    Returns
    -------
    stdout : str
    stderr : str
    exit_code : int
    """
    pp = _sub.Popen(command, shell=True, universal_newlines=True,
                    stdout=_sub.PIPE, stderr=_sub.PIPE)
    out, err = pp.communicate()
    code = pp.returncode
    if raise_on_error and code != 0:
        raise _sub.CalledProcessError(
            returncode=code, cmd=command, output=out, stderr=err
        )
    return out.decode(), err.decode(), code


class PLINK(object):

    """A reusable object to run plink jobs quickly."""

    written_files = {}

    def __init__(self, pop_file=None):
        """Load population information."""
        if not pop_file:
            pop_file = _os.path.join(
                DATA_DIR, 'integrated_call_samples_v3.20130502.ALL.panel'
            )
        individuals = {}
        with open(pop_file) as fin:
            assert fin.readline() == 'sample\tpop\tsuper_pop\tgender\t\t\n'
            for line in fin:
                ind, pop, _, _ = line.split('\t')
                if pop not in individuals:
                    individuals[pop] = []
                individuals[pop].append(ind)
        self.individuals = individuals

    def pop_file(self, populations=None):
        """Write temp file with a list of individuals in population."""
        populations = populations if populations else individuals.keys()
        if isinstance(populations, str):
            populations = [populations]
        populations = list(populations)
        if ','.join(populations) in self.written_files:
            fl = self.written_files[populations]
            if _os.path.isfile(fl):
                return fl
        pop_ids  = []
        bad_pops = []
        for pop_i in populations:
            if pop_i in individuals:
                pop_ids += individuals[pop_i]
            else:
                bad_pops.append(pop_i)
        if bad_pops:
            err = (
                "{} are not ancestral populations. Choose one of the following "
                "ancestral populations: {}"
            ).format(bad_pops, POPULATIONS)
            raise ValueError(err)

        pop_ids = sorted(set(pop_ids))

        _, file_location = _temp(prefix='-'.join(populations), dir='/tmp')
        with open(file_location, 'w') as outfile:
            outfile.write('\n'.join(pop_ids))

        self.written_files[','.join(populations)] = file_location
        return file_location

    def one_to_many(self, snp, comp_list, chrom, r2=0.9, populations=None):
        """Get one-to-many LD information using plink.

        Parameters
        ----------
        snp : str
            rsID of a SNP to query
        comp_list : list_of_str
            list of rsIDs to compare to
        chrom : str
            which chromosome to search
        r2 : float
            r-squared level to use for filtering
        populations : list_of_str
            list of populations to include in the analysis

        Returns
        -------
        matching_snps
            For every matching SNP that beats the r-squared: {
                snp: {r2: r-squared, dprime: d-prime, phased: phased-alleles}
            }
        """
        _, temp_file = _temp(prefix='plink', dir='/tmp')
        pop_file = self.pop_file(populations)
        # We need to include the query SNP in the lookup list
        comp_list.append(snp)
        comp_list = sorted(set(comp_list))
        # Build the command
        plink_cmnd = (
            '{plink} --bfile {bfile} --r2 in-phase dprime --ld-snp {snp} '
            '--snps {comp_list} --keep {ind_file} --out {tmp}'
        ).format(
            plink='plink',
            bfile=_os.path.join(
                DATA_DIR,
                'ALL.{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes'.format(chrom)
            ),
            snp=snp, comp_list=' '.join(comp_list), pop_file=pop_file,
            tmp=temp_file
        )
        # Run it
        run(plink_cmnd, True)
        # Parse the output file
        results = {}
        with open(temp_file + '.ld') as fin:
            for line in fin:
                f = _re.split(r' +', line.strip())
                snp2, phased, rsquared, dprime = f[5], f[6], f[7], f[8]
                if snp2 == snp:
                    continue
                results[snp2] = {'r2': rsquared, 'dprime': dprime,
                                    'phased': phased}
        _os.remove(temp_file + '.ld')
        return results


# Primary function
def comp_snp_lists(list1, list2, populations=None, r2=0.9, distance='50kb',
                   return_dataframe=False, raise_on_error=False):
    """Compare two SNP lists, filter by r2.

    Runs a dbSNP lookup to get rsID and location information for all SNPs in
    list1 and list2, then creates pairs of SNPs: for every SNP in list1, create
    a pair with every SNP in list2 within 'distance'. Then, do a one-to-many
    LD comparison with plink for every snp-list combination. This calculation
    is filtered by population.

    Parameters
    ----------
    list1/list2 : list_of_str
        list of rsIDs and chr:pos (where pos is base-1)
    populations : list_of_str
        list of 1000genomes populations, if not provided, all possible
        populations are used. This list is used for the LD calculations and
        limits the individuals considered in the LD calculation.
    r2 : float
        An r-squared cutoff to use for filtering
    distance : str_or_int
        A distance away from a SNP in list1 to consider SNPs in list2, can be
        a integer of bases or a string with a suffix kb (kilo-bases) or mb
        (mega-bases)
    return_dataframe : bool
        Return a dataframe instead of a list of SNP_Link objects
    raise_on_error : bool, optional
        Raise an Exception if there is a problem with any of the SNPs in the
        list, otherwise errors are silently ignored.

    Returns
    -------
    dict_or_DataFrame
        If return_dataframe is True, then a DataFrame is built and returned
        otherwise a dictionary of {SNP: [SNP_Link, ...]}.
    """
    snps1 = list_to_rsid_and_locs(list1, raise_on_error)
    snps2 = list_to_rsid_and_locs(list2, raise_on_error)
    pairs = create_pairs(snps1, snps2)
    pairs = filter_pairs(pairs, r2, populations)


def filter_pairs(pairs, r2=0.9, populations=None):
    """Use plink to lookup LD and filter lists by the r2."""
    plink = PLINK()
    run_lists = {}
    for chrom, data in pairs.items():
        run_lists[chrom] = {}
        for snp, info in data.items():
            if info['matches']:
                run_lists[chrom][snp] = info['matches']
    results = {}
    for chrom, snps in run_lists.items():
        for snp, comp_list in snps.items():
            results[snp] = plink.one_to_many(
                snp, comp_list, chrom, r2, populations
            )
    # Add all data back again
    filtered = {}
    for chrom, data in pairs.items():
        for snp, info in data.items():
            filtered[snp] = {'loc': info['loc']}
            if info['matches']:
                filtered[snp]['matches'] = {}
                for snp2, pos2 in info['matches']:
                    if snp2 in results[snp]:
                        filtered[snp]['matches'][snp2] = results[snp][snp2]
                        filtered[snp]['matches'][snp2]['loc'] = pos2
    return filtered


def create_pairs(list1, list2, distance='50kb'):
    """Create a dictionary of pairwise combinations from list1 and list2.

    For every SNP in list1, checks all SNPs in list2 to see which are within
    distance of SNP1.

    Parameters
    ----------
    list1/list2 : list_of_tuples
        lists from list_to_rsid_and_locs() in the format:
        [(rsid, chromosome, int(position))]  # position is base-1
    distance : str_or_int
        A distance away from a SNP in list1 to consider SNPs in list2, can be
        a integer of bases or a string with a suffix kb (kilo-bases) or mb
        (mega-bases)

    Returns
    -------
    dict
        A dictionary of {chrom:
            SNP1: {
                loc: position,
                matches: [(snp2, position2), ...]
            }
        }
    """
    if isinstance(distance, str):
        dist, mod = _re.match(r'(\d+)(\D+)', '50kb').groups()
        dist = int(dist)
        mod = mod.lower()
        if not mod in ['kb', 'mb']:
            raise ValueError('Cannot parse {}, must be in kb or mb only'
                             .format(distance))
        if mod == 'kb':
            dist = dist * 1000
        elif mod == 'mb':
            dist = dist * 1000000
    # Covert list2 to a dictionary for speed
    dict2 = {}
    for name, chrom, loc in list2:
        if chrom not in dict2:
            dict2[chrom] = []
        dict2[chrom].append((loc, name))
    results = {}
    for snp1, chrom1, pos1 in list1:
        if chrom1 not in results:
            results[chrom1] = {}
        results[chrom1][snp1] = {'loc': pos1, 'matches': []}
        if chrom1 not in dict2:
            continue
        for pos2, snp2 in dict2[chrom1]:
            if snp1 == snp2:
                continue
            if abs(pos2 - pos1) > dist:
                continue
            results[chrom1][snp1]['matches'].append((snp2, pos2))
    return results


def list_to_rsid_and_locs(list1, raise_on_missing=False):
    """Convert every entry in each list to an rsid and chr:loc position.

    dbSNP is base-0 but the chr:loc positions are all treated as base-1, so any
    chr:loc positions should be base-1.

    Uses a single large dbSNP lookup to create a dictionary and then parses
    the two lists.

    Parameters
    ----------
    list1 : list_of_str
        list of rsIDs and chr:pos (where pos is base-1)
    raise_on_missing : str
        Raise an Exception if there is a problem with any of the SNPs (i.e.
        they aren't in dbSNP or are indels

    Raises
    ------
    MissingSNPError
        If a SNP is missing from dbSNP
    BadSNPError
        If a SNP is an indel

    Returns
    -------
    list_of_tuples
        [(rsid, chromosome, int(position))]  # position is base-1
    """
    list1_rs = []
    list1_locs = []
    bad = []
    for snp in list1:
        if snp.startswith('rs'):
            list1_rs.append(snp)
        elif ':' in snp:
            list1_locs.append(snp)
        else:
            bad.append(snp)
    if bad:
        err = (
            'The following SNPs fail format requirements:\n{}'
            .format(failed)
        )
        if raise_on_missing:
            raise MissingSNPError(err)
        _sys.stderr.write(err + '\n')

    if list1_rs:
        db = _dbSNP.DB(DB_PATH, DB_VERS)
        rs_lookup = {
            i[0]: (i[0], i[1], i[2]+1) for i in db.query(
                db.Row.name, db.Row.chrom, db.Row.start, db.Row.end
            ).filter(
                db.Row.name.in_(list1_rs)
            )
            if not i[3]-i[2] > 1
        }

    if list1_locs:
        db = _dbSNP.DB(DB_PATH, DB_VERS)
        query = {}
        for chrom, loc in list1_locs.split(':'):
            if chrom not in query:
                query[chrom] = []
            query[chrom].append(int(loc-1))
        loc_lookup = {
            '{}:{}'.format(i.chrom, i.start+1): (i.name, i.chrom, i.start)
            for i in db.lookup_locations(query)
            if not i.end-i.start > 1
        }

    failed = []
    new_lst = []
    for snp in list1:
        if snp.startswith('rs'):
            if not snp in rs_lookup:
                err = '{} not in dbSNP'.format(snp)
                if raise_on_missing:
                    raise MissingSNPError(err)
                else:
                    sys.stderr.write(err + '\n')
            new_lst.append(rs_lookup[snp])
        elif ':' in snp:
            if not snp in loc_lookup:
                err = '{} not in dbSNP'.format(snp)
                if raise_on_missing:
                    raise MissingSNPError(err)
                else:
                    sys.stderr.write(err + '\n')
            new_lst.append(loc_lookup[snp])
    return new_lst

def get_arg_parser():
    """Create an argument parser."""
    parser  = _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter)

    # Positional arguments
    parser.add_argument(
        'list1',
        help="SNP list 1, newline separated list of rsids or chr:loc"
    )
    parser.add_argument(
        'list2',
        help="SNP list 2, newline separated list of rsids or chr:loc"
    )
    parser.add_argument(
        'outfile', nargs='?',
        help="File to write output to, default is STDOUT"
    )

    # Filtering options
    filtering = parser.add_argument_group('filter', 'Filtration options')
    filtering.add_argument('--r2', default=0.9,
                           help='Minimum r-squared to consider LD (0.9)')
    filtering.add_argument('--distance', default='50kb',
                           help='Max distance to consider LD (50kb)')

    # Optional flags
    parser.add_argument('-p', '--pandas', action="store_true",
                        help="Write file as a pandas DataFrame instead.")
    parser.add_argument('--populations',
                        help="Comma separated list of populations to check " +
                        "Default: all populations.")

    return parser


def main(argv=None):
    """Run as a script."""
    if not argv:
        argv = _sys.argv[1:]

    parser = get_arg_parser()

    args = parser.parse_args(argv)

    if not args.outfile and args.pandas:
        __sys.stderr.write('Cannot write pandas DataFrame to STDOUT\n')
        return 1

    with open(args.list1) as fin:
        list1 = fin.read().strip().split('\n')
    with open(args.list2) as fin:
        list2 = fin.read().strip().split('\n')

    if args.populations:
        pops = args.populations.split(',')
    else:
        pops = None

    if args.pandas:
        table_args = {'output_table': None, 'return_dataframe': True}
    else:
        table_args = {'output_table': args.outfile, 'return_dataframe': False}

    out = comp_snp_lists(list1, list2, populations, r2=float(args.r2),
                         distance=args.distance, **table_args)

    if args.pandas:
        out.to_pickle(args.outfile)

if __name__ == '__main__' and '__file__' in globals():
    _sys.exit(main())
