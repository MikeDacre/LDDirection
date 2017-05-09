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
import pickle as _pickle
import subprocess as _sub
import argparse as _argparse
import multiprocessing as _mp
import bz2 as _bz2
import gzip as _gzip
from tempfile import mkstemp as _temp

try:
    import pandas as _pd
except ImportError:
    _pd = None

try:
    import fyrd as _fyrd
except ImportError:
    _fyrd = None

import dbSNP as _dbSNP

from .snp_link import SNP_Pair

try:
    from tqdm import tqdm, tqdm_notebook
    try:
        if str(type(get_ipython())) == "<class 'ipykernel.zmqshell.ZMQInteractiveShell'>":
            pb = tqdm_notebook
        else:
            pb = tqdm
    except NameError:
        pb = tqdm
except ImportError:
    pb = None

from .distance_filtering import distance_filter

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
        _sys.stderr.write(
            '{}\nfailed with:\nCODE: {}\nSTDOUT:\n{}\nSTDERR:\n{}\n'
            .format(command, code, out, err)
        )
        raise _sub.CalledProcessError(
            returncode=code, cmd=command, output=out, stderr=err
        )
    return out, err, code


def _chrom_sort(x):
    """Return an integer for sorting from a chromosome."""
    if x.startswith('chr'):
        x = x[3:]
    if x.upper() == 'X':
        return 100
    elif x.upper() == 'Y':
        return 101
    elif x.upper().startswith('M'):
        return 150
    elif x.isdigit():
        return int(x)
    else:
        return x


class PLINK(object):

    """A reusable object to run plink jobs quickly."""

    written_files = {}
    bims = {}

    def __init__(self, pop_file=None, plink='plink'):
        """Load population information.

        Parameters
        ----------
        pop_file : str, optional
            Path to the 1000 genomes population file
        plink : str, optional
            Path to plink executable, otherwise searches PATH
        """
        self.plink = plink
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
        populations = populations if populations else self.individuals.keys()
        if isinstance(populations, str):
            populations = [populations]
        populations = list(populations)
        if ','.join(populations) in self.written_files:
            fl = self.written_files[','.join(populations)]
            if _os.path.isfile(fl):
                return fl
        pop_ids  = []
        bad_pops = []
        for pop_i in populations:
            if pop_i in self.individuals:
                pop_ids += self.individuals[pop_i]
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
            outfile.write('\n'.join(
                [' '.join(x) for x in zip(pop_ids, pop_ids)]
            ))

        self.written_files[','.join(populations)] = file_location
        return file_location

    def bim_snps(self, chrom):
        """Return and cache all SNPs in a bim file."""
        if chrom in self.bims:
            with open(self.bims[chrom], 'rb') as fin:
                return _pickle.load(fin)
        bfile = _os.path.join(
            DATA_DIR,
            'ALL.{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes'
            .format(chrom)
        )
        bim = bfile + '.bim'
        rsids = []
        with open(bim) as fin:
            for line in fin:
                rsids.append(line.split('\t')[1])
        rsids = frozenset(rsids)
        _, pfl = _temp()
        with open(pfl, 'wb') as fout:
            _pickle.dump(rsids, fout)
        self.bims[chrom] = pfl
        return rsids

    def one_to_many(self, snp, comp_list, chrom, r2=0.6, populations=None,
                    raise_on_error=False, logfile=_sys.stderr):
        """Get one-to-many LD information using plink.

        Parameters
        ----------
        snp : str
            rsID of a SNP to query
        comp_list : list_of_str
            list of rsIDs to compare to
        chrom : str
            which chromosome to search
        r2 : float, optional
            r-squared level to use for filtering
        populations : list_of_str, optional
            list of populations to include in the analysis
        raise_on_error : bool, optional
            if False, will return None if primary SNP missing from bim file
        logfile : filehandle, optional
            A file like object to write to

        Returns
        -------
        matching_snps : dict
            For every matching SNP that beats the r-squared: {
                snp: {r2: r-squared, dprime: d-prime, phased: phased-alleles}
            }
            If plink job fails, returns an empty dictionary.
        """
        _, temp_file = _temp(prefix='plink', dir='/tmp')
        pop_file = self.pop_file(populations)
        bfile = _os.path.join(
            DATA_DIR,
            'ALL.{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes'
            .format(chrom)
        )
        bim = bfile + '.bim'
        # We need to include the query SNP in the lookup list
        comp_list = list(comp_list)
        comp_list.append(snp)
        comp_list = sorted(set(comp_list))
        # Filter SNPs not in the bim
        bim_snps = self.bim_snps(chrom)
        bad = []
        if snp not in bim_snps:
            err = ('Primary SNP {} is not in BIM {}, cannot continue.'
                   .format(snp, bim))
            if raise_on_error:
                raise BadSNPError(err)
            else:
                logfile.write(err + '\n')
                return None
        bad = []
        for s in comp_list:
            if s not in bim_snps:
                bad.append(s)
                comp_list.remove(s)
        if bad:
            _sys.stderr.write(('{} removed from comparison list as not in ' +
                               'bim file\n').format(bad))
        del(bim_snps)
        # Build the command
        plink_cmnd = (
            '{plink} --bfile {bfile} --r2 in-phase dprime --ld-snp {snp} '
            '--snps {comp_list} --keep {ind_file} --out {tmp}'
        ).format(
            plink=self.plink,
            bfile=bfile,
            snp=snp, comp_list=' '.join(comp_list),
            ind_file=pop_file, tmp=temp_file
        )
        # Run it
        stdout, stderr, code = run(plink_cmnd, raise_on_error)
        # Parse the output file
        if code != 0:
            logfile.write(
                '{}: plink command failed'.format(snp) +
                'Command: {}\nExit Code: {}\nSTDOUT:\n{}\bSTDERR:\n{}\n'
                .format(plink_cmnd, code, stdout, stderr)
            )
            return {}
        results = {}
        with open(temp_file + '.ld') as fin:
            # Check header
            line = fin.readline().strip()
            assert _re.split(r' +', line) == [
                'CHR_A', 'BP_A', 'SNP_A', 'CHR_B', 'BP_B',
                'SNP_B', 'PHASE', 'R2', 'DP'
            ]
            for line in fin:
                f = _re.split(r' +', line.strip())
                snp2, phased = f[5], f[6]
                rsquared, dprime = float(f[7]), float(f[8])
                if snp2 == snp:
                    continue
                if rsquared < r2:
                    continue
                try:
                    p1, p2 = phased.split('/')
                    s1a, s2a = p1
                    s1b, s2b = p2
                    lookup = {snp:  {s1a: s2a, s1b: s2b},
                              snp2: {s2a: s1a, s2b: s1b}}
                except ValueError:
                    lookup = {}
                results[snp2] = {'r2': rsquared, 'dprime': dprime,
                                 'phased': phased, 'lookup': lookup}
        _os.remove(temp_file + '.ld')
        return results


###############################################################################
#                              dbSNP Connection                               #
###############################################################################


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
        Note, itmes can also be a three-tuple:
            (rsid, chrom, loc)
        In this format, no lookup is performed, so function runs very fast.
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
    dict_of_lists
        {'chrom': [(rsid, int(position))]}  # position is base-1
    """
    list1_rs = []
    list1_locs = []
    list1_done = {}
    bad = []
    for snp in list1:
        if isinstance(snp, str) and '\t' in snp:
            snp = snp.strip().split('\t')
        if isinstance(snp, (list, tuple)):
            if not len(snp) == 3:
                bad.append(snp)
            else:
                rsid, chrom, loc = snp
                chrom = chrom if chrom.startswith('chr') else 'chr' + chrom
                loc   = int(loc)
                list1_done[rsid] = (chrom, loc)
        elif snp.startswith('rs'):
            list1_rs.append(snp)
        elif ':' in snp:
            list1_locs.append(snp)
        else:
            bad.append(snp)
    if bad:
        err = (
            'The following SNPs fail format requirements:\n{}'
            .format(bad)
        )
        if raise_on_missing:
            raise MissingSNPError(err)
        _sys.stderr.write(err + '\n')

    if list1_rs:
        _sys.stderr.write('Doing rsID DB lookup.\n')
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
        _sys.stderr.write('Doing position DB lookup.\n')
        db = _dbSNP.DB(DB_PATH, DB_VERS)
        query = {}
        for chrom, loc in [i.split(':') for i in list1_locs]:
            if chrom not in query:
                query[chrom] = []
            query[chrom].append(int(loc-1))
        loc_lookup = {
            '{}:{}'.format(i.chrom, i.start+1): (i.name, i.chrom, i.start+1)
            for i in db.lookup_locations(query)
            if not i.end-i.start > 1
        }

    failed  = []
    results = {}
    done = []
    for snp in list1:
        if isinstance(snp, str) and '\t' in snp:
            snp = snp.strip().split('\t')
        if isinstance(snp, (list, tuple)):
            if not snp[0] in list1_done:
                failed.append(snp[0])
            s_rs = snp[0]
            s_chrom, s_loc = list1_done[s_rs]
        elif snp.startswith('rs'):
            if not snp in rs_lookup:
                failed.append(snp)
                continue
            s_rs, s_chrom, s_loc = rs_lookup[snp]
        elif ':' in snp:
            if not snp in loc_lookup:
                failed.append(snp)
                continue
            s_rs, s_chrom, s_loc = loc_lookup[snp]
        if s_chrom not in results:
            results[s_chrom] = []
        if s_rs in done:
            continue
        done.append(s_rs)
        results[s_chrom].append((s_rs, s_loc))
    if failed:
        err = '{} not in dbSNP'.format(snp)
        if raise_on_missing:
            raise MissingSNPError(err)
        else:
            _sys.stderr.write(err + '\n')
    return results


def save_list(snp_list, outfile=None, pickle=False):
    """Convert a list of rsIDs or positions to a pre-processed format.

    Converts either rsID or chr:loc (in base-1) style SNP locations into
    (rsid, chrom, loc) (also base-1) to make future lookups faster.

    Parameters
    ----------
    snp_list : list_of_str
        list of rsIDs and chr:pos (where pos is base-1)
        Note, itmes can also be a three-tuple:
            (rsid, chrom, loc)
    outfile : str, optional
        A file to write the output to.
    pickle : bool, optional
        Write the tuple list as a pickled object instead of text.

    Returns
    -------
    list
        A list of (rsid, chrom, loc)

    Writes
    ------
    rsid\tchrom\tloc
    """
    _sys.stderr.write('Parsing input list.\n')
    parsed_list = list_to_rsid_and_locs(snp_list)
    _sys.stderr.write('Done.\n')
    output = []
    for chrom, info in parsed_list.items():
        for rsid, loc in info:
            output.append((rsid, chrom, loc))
    if outfile:
        _sys.stderr.write('Writing output.\n')
        if pickle:
            with open_zipped(outfile, 'wb') as fout:
                _pickle.dump(output, fout)
        else:
            with open_zipped(outfile, 'w') as fout:
                for row in output:
                    rsid, chrom, loc = row
                    fout.write('{}\t{}\t{}\n'.format(rsid, chrom, loc))
    return output


def open_zipped(infile, mode='r'):
    """Return file handle of file regardless of compressed or not.

    Also returns already opened files unchanged, text mode automatic for
    compatibility with python2.
    """
    # Return already open files
    if hasattr(infile, 'write'):
        return infile
    # Make text mode automatic
    if len(mode) == 1:
        mode = mode + 't'
    if not isinstance(infile, str):
        raise ValueError("I cannot open a filename that isn't a string.")
    if infile.endswith('.gz'):
        return _gzip.open(infile, mode)
    if infile.endswith('.bz2'):
        if hasattr(_bz2, 'open'):
            return _bz2.open(infile, mode)
        else:
            return _bz2.BZ2File(infile, mode)
    return open(infile, mode)


###############################################################################
#                             Filtering Functions                             #
###############################################################################


def filter_by_distance(ref_snps, comp_snps, distance='50kb'):
    """Create a dictionary of pairwise combinations from list1 and list2.

    For every SNP in ref_snps, checks all SNPs in comp_snps to see which are
    within distance of SNP1.

    Parameters
    ----------
    ref_snps/comp_snps : dict_of_lists
        In the format: {chrom: [(rsid, location), ...]} where location is
        base-1. This format is output by list_to_rsid_and_locs().
    distance : str_or_int, optional
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
    # Calculate distance
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

    results = {}
    for chrom, snps in ref_snps.items():
        results[chrom] = distance_filter(snps, comp_snps[chrom], dist)

    return results


def _ld_job(snps, plink, chrom, r2, populations, pbar=True):
    """Private function to run plink, see filter_by_ld()."""
    if pbar and pb:
        it = pb(snps, unit='plink_calculations')
        log = it
    else:
        it = snps
        log = _sys.stderr
    snp_results = {}
    for snp, _, comp_list in it:
        if not comp_list:
            continue
        comp_snps = [rsid for rsid, loc in comp_list]
        snp_results[snp] = plink.one_to_many(
            snp, comp_snps, chrom, r2, populations, raise_on_error=False,
            logfile=log
        )
    return snp_results


def filter_by_ld(pairs, r2=0.6, populations=None, plink='plink'):
    """Use plink to lookup LD and filter lists by the r2.

    Parameters
    ----------
    r2 : float
        An r-squared cutoff to use for filtering
    distance : str_or_int
        A distance away from a SNP in list1 to consider SNPs in list2, can be
        a integer of bases or a string with a suffix kb (kilo-bases) or mb
        (mega-bases)
    plink : str, optional
        Path to plink executable, otherwise searches PATH

    """
    plink = PLINK(plink=plink)
    l = 0
    for v in pairs.values():
        for i in [i[2] for i in v]:
            l += len(i)
    multi = True if l > 200 else False
    if multi:
        jobs = {}
        if _fyrd and _fyrd.queue.MODE != 'local':
            print('Running jobs on the cluster with fyrd.')
        else:
            pool = _mp.Pool()
            print('Running jobs in parallel')

    results = {}
    for chrom in sorted(pairs, key=_chrom_sort):
        snps = pairs[chrom]
        print('Working on chromosome {}'.format(chrom))
        args = (snps, plink, chrom, r2, populations)
        if multi:
            args += (False, )
            if _fyrd and _fyrd.queue.MODE != 'local':
                jobs[chrom] = _fyrd.submit(
                    _ld_job, args, walltime="01:00:00", mem='8G', cores=4)
            else:
                jobs[chrom] = pool.apply_async(
                    _ld_job, args)
        else:
            results[chrom] = _ld_job(*args)

    if multi:
        print('Getting results from parallel jobs')
        for chrom, job in jobs.items():
            results[chrom] = job.get()
        if not _fyrd and _fyrd.queue.MODE != 'local':
            pool.close()

    # Add all data back again
    filtered = {}
    for chrom, data in pairs.items():
        for snp, loc, comp_list in data:
            filtered[snp] = {'loc': loc, 'chrom': chrom}
            if snp in results[chrom]:
                if results[chrom][snp]:
                    filtered[snp]['matches'] = {}
                    for snp2, pos2 in comp_list:
                        if snp2 in results[chrom][snp]:
                            filtered[snp]['matches'][snp2] = results[chrom][snp][snp2]
                            filtered[snp]['matches'][snp2]['loc'] = pos2
    return filtered


###############################################################################
#                          Output Parsing Functions                           #
###############################################################################


def pairs_to_dataframe(pairs):
    """Convert output of filter_by_ld() to a dataframe."""
    if not _pd:
        raise ImportError('pandas is not available')
    header = ['SNP', 'Chrom', 'Location', 'Match', 'Match_Loc',
              'SNP1_A1', 'SNP1_A2', 'SNP2_A1', 'SNP2_A2',
              'R2', 'Dprime', 'SNP1_A1_Pair']
    rows  = []
    index = []
    for snp1, data in pairs.items():
        cols = [snp1, data['chrom'], data['loc']]
        if not 'matches' in data:
            continue
        for snp2, info in data['matches'].items():
            s1_alleles = sorted(info['lookup'][snp1])
            s2_alleles = [info['lookup'][snp1][i] for i in s1_alleles]
            rcols = [
                snp2, info['loc'],
                s1_alleles[0], s1_alleles[1], s2_alleles[0], s2_alleles[1],
                info['r2'], info['dprime'], s2_alleles[0]
            ]
            rows.append(cols + rcols)
            index.append('{}:{}'.format(snp1, snp2))
    return _pd.DataFrame(rows, columns=header, index=index).sort_values(
        ['Chrom', 'Location'])


def pairs_to_SNP_Pairs(pairs, populations):
    """Convert output of filter_by_ld() to dict of SNP_Pairs."""
    results = {}
    for snp1, info in pairs.items():
        results[snp1] = []
        if info['matches']:
            for snp2, data in info['matches'].items():
                results[snp1].append(
                    SNP_Pair(
                        plink = {
                            snp1: {
                                'chrom': info['chrom'],
                                'loc' : info['loc'],
                                'matches': {
                                    snp2: data
                                }
                            }
                        },
                        populations = populations
                    )
                )
    return results


###############################################################################
#                     Primary Function to Run Everything                      #
###############################################################################


def comp_snp_lists(list1, list2, populations=None, r2=0.6, distance='50kb',
                   plink='plink', return_dataframe=False, return_snppair=False,
                   raise_on_error=False):
    """Compare two SNP lists, filter by r2.

    Runs a dbSNP lookup to get rsID and location information for all SNPs in
    list1 and list2, then creates pairs of SNPs: for every SNP in list1, create
    a pair with every SNP in list2 within 'distance'. Then, do a one-to-many
    LD comparison with plink for every snp-list combination. This calculation
    is filtered by population.

    Parameters
    ----------
    list1/list2 : list_of_str or list_of_tuple
        list of rsIDs and chr:pos (where pos is base-1)
        Can also be a list of three-tuples in the format:
            (str(rsid), str(chrom), int(loc))
        In this form, no database lookup is performed, making the query much
        faster and removing the need to the SNP to be in dbSNP.
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
    plink : str, optional
        Path to plink executable, otherwise searches PATH
    return_dataframe : bool, optional
        Return a dataframe instead of JSON
    return_snppair : bool, optional
        Return a dictionary of SNP_Link objects
    raise_on_error : bool, optional
        Raise an Exception if there is a problem with any of the SNPs in the
        list, otherwise errors are silently ignored.

    Returns
    -------
    JSON : dict
        If return_dataframe and return_snppairs are False, a json-style
        dictionary is returned with the format:
            {snp1: {chrom: chrX, loc: #, matches:
                {snp2: {loc: #, r2: float, dprime: float, phased: AT/GC,
                        lookup: {snp1: {S1A1: S2A1, ..}, snp2: {..}}}}
            }}
    DataFrame : pandas.DataFrame
        If return_dataframe is True, then a DataFrame is built and returned
        with the rows:
            SNP, Chrom, Location, Match, Match_Loc, SNP1_A1, SNP1_A2,
            SNP2_A1, SNP2_A2, R2, Dprime, SNP1_A1_Pair
        The SNPs are in order, so SNP1_A1 occurs with SNP2_A2, etc.
    SNP_Pairs : dict_of_SNP_Pairs
        If return_snplink is True, a dictionary of SNP_Link objects is returned
            {snp1: [SNP_Pair, SNP_Pair, ...]}
    """
    print('Getting SNP info for list1')
    snps1 = list_to_rsid_and_locs(list1, raise_on_error)
    print('Getting SNP info for list2')
    snps2 = list_to_rsid_and_locs(list2, raise_on_error)
    print('Getting all possible pairs within {}'.format(distance))
    pairs = filter_by_distance(snps1, snps2, distance=distance)
    l = 0
    for v in pairs.values():
        for i in [i[2] for i in v]:
            l += len(i)
    print('Running plink filters on {} calculations'.format(l))
    pairs = filter_by_ld(pairs, r2, populations=populations, plink=plink)
    print('Done')
    if return_dataframe:
        return pairs_to_dataframe(pairs)
    if return_snppair:
        return pairs_to_SNP_Pairs(pairs, populations=populations)
    return pairs


###############################################################################
#                               Run as a Script                               #
###############################################################################

analyze_usage = """\
ldlists analyze --r2 0.8 --distance 100kb snp_list comp_list [outfile]
"""

pre_usage = """\
ldlists preprocess snp_list [outfile]
"""

general_usage = """
General Usage
-------------
ldlists [--help] <mode> <options> snp_list [comp_list] [outfile]

Examples
--------
ldlists --help  # Prints detail about algorithm
{}
{}

""".format(analyze_usage.strip(), pre_usage.strip())

analyze_desc = """\
Run the full analysis on two lists of SNPs.
"""

analyze_epilog = """\
snp_list and comp_list must be files containing newline separated SNP info.

Each line can be one of the following:
    - an rsID
    - a chr:loc style name, must be base-1
    - a tab delimited line with three columns:
        rsid, chromosome, location (base-1)
If the line is a tab-delimited line, no database lookup is done, saving a lot
of time for large lists. Otherwise a dbSNP lookup is required.

To create a file in this format, run this program in preprocess mode.
"""

pre_desc = """\
Convert a list of SNP info into a complete set of information to allow running
in analyze mode without doing a DB lookup.
"""

pre_epi = """\
Input file must be a newline separated file where each line can be one of:
    - an rsID
    - a chr:loc style name, must be base-1
    - a tab delimited line with three columns:
        rsid, chromosome, location (base-1)

The output file will be three column tab-delimited file of:
    rsid, chromosome, location (base-1)
"""

def get_arg_parser():
    """Create an argument parser."""
    parser  = _argparse.ArgumentParser(
        description=__doc__, #usage=general_usage,
        formatter_class=_argparse.RawDescriptionHelpFormatter)

    subparsers = parser.add_subparsers(title='mode', dest='mode')

    # Primary mode
    analyze = subparsers.add_parser(
        'analyze', usage=analyze_usage,
        description=analyze_desc, epilog=analyze_epilog,
        formatter_class=_argparse.RawDescriptionHelpFormatter)

    # Positional arguments
    analyze.add_argument(
        'snp_list',
        help="File containing primary list of SNPs to analyze"
    )
    analyze.add_argument(
        'comp_list',
        help="File containing list of SNPs to compare to"
    )
    analyze.add_argument(
        'outfile', nargs='?',
        help="File to write output to, default is STDOUT"
    )

    # Filtering options
    filtering = analyze.add_argument_group('filter', 'Filtration options')
    filtering.add_argument('--r2', default=0.9,
                           help='Minimum r-squared to consider LD (0.9)')
    filtering.add_argument('--distance', default='50kb',
                           help='Max distance to consider LD (50kb)')

    # Optional flags
    analyze.add_argument('-p', '--pandas', action="store_true",
                         help="Write file as a pandas DataFrame instead.")
    analyze.add_argument('--populations',
                         help="Comma separated list of populations to check " +
                         "Default: all populations.")

    # Dump Lists
    pre = subparsers.add_parser(
        'preprocess', usage=pre_usage, description=pre_desc, epilog=pre_epi,
        formatter_class=_argparse.RawDescriptionHelpFormatter)

    pre.add_argument(
        'snp_list',
        help="File containing primary list of SNPs to parse"
    )
    pre.add_argument(
        'outfile', nargs='?',
        help="File to write output to, default is STDOUT"
    )

    return parser


def main(argv=None):
    """Run as a script."""
    if not argv:
        argv = _sys.argv[1:]

    parser = get_arg_parser()

    args = parser.parse_args(argv)

    if not args.mode:
        parser.print_usage()
        _sys.stdout.write(general_usage)
        _sys.stderr.write('Mode required.\n')
        return 2

    with open_zipped(args.snp_list) as fin:
        snp_list = fin.read().strip().split('\n')

    if args.mode == 'preprocess':
        outfile = args.outfile if args.outfile else _sys.stdout
        save_list(snp_list, outfile)
        return 0

    if not args.outfile and args.pandas:
        _sys.stderr.write('Cannot write pandas DataFrame to STDOUT\n')
        return 1

    with open_zipped(args.comp_list) as fin:
        comp_list = fin.read().strip().split('\n')

    if args.populations:
        pops = args.populations.split(',')
    else:
        pops = None

    out = comp_snp_lists(snp_list, comp_list, pops, r2=float(args.r2),
                         distance=args.distance, return_dataframe=True)

    if args.pandas:
        out.to_pickle(args.outfile)
    else:
        out.to_csv(args.outfile, sep='\t')

    return 0

if __name__ == '__main__' and '__file__' in globals():
    _sys.exit(main())
