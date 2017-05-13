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
import sqlite3 as _sqlite3
import pickle as _pickle
import subprocess as _sub
import argparse as _argparse
import multiprocessing as _mp
import bz2 as _bz2
import gzip as _gzip
from time import sleep as _sleep
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

MIN_BUNDLE = 150  # Minimum number of jobs per node, 150 takes less than a minute.

DB_PATH = '/godot/dbsnp/db/hg19'
DB_VERS = 150

POPULATIONS = ["ALL", "AFR", "AMR", "EAS", "EUR", "SAS", "ACB", "ASW", "BEB",
               "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH",
               "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL",
               "PJL", "PUR", "STU", "TSI", "YRI"]

# Set data directories
DATA_DIR = "/godot/1000genomes/1000GP_Phase3"

__all__ = ['comp_snp_lists', 'list_to_rsid_and_locs', 'filter_by_distance',
           'filter_by_ld', 'save_list']


###############################################################################
#                   Run plink jobs without repetative code                    #
###############################################################################


class PLINK(object):

    """A reusable object to run plink jobs on 1000genomes files quickly."""

    written_files = {}
    bims = {}

    def __init__(self, data=None, pop_file=None, plink='plink'):
        """Load population information.

        Parameters
        ----------
        data : str, optional
            Path the the 1000genomes data files
            (e.g. ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.{bed,bim})
            Default set in DATA_DIR hardcoded in this file.
        pop_file : str, optional
            Path to the 1000 genomes population file
        plink : str, optional
            Path to plink executable, otherwise searches PATH
        """
        self.plink = plink
        data = data if data else DATA_DIR
        if not _os.path.isdir(data):
            raise ValueError('{} must be a directory, but no directory found'
                             .format(data))
        self.path = _os.path.abspath(data)
        if not pop_file:
            pop_file = _os.path.join(
                self.path, 'integrated_call_samples_v3.20130502.ALL.panel'
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
        self.files = {}
        for fl in _os.listdir(self.path):
            if 'phase3_' not in fl:
                continue
            if not fl.endswith('bed'):
                continue
            if not fl.startswith('ALL'):
                continue
            chrom = fl.split('.')[1]
            if not chrom.startswith('chr'):
                continue
            fl    = _os.path.abspath(_os.path.join(self.path, fl))
            root  = '.'.join(fl.split('.')[:-1])
            bed   = _os.path.abspath('{}.bed'.format(root))
            bim   = _os.path.abspath('{}.bim'.format(root))
            fam   = _os.path.abspath('{}.fam'.format(root))
            assert _os.path.isfile(bed)
            assert _os.path.isfile(bim)
            assert _os.path.isfile(fam)
            c1 = chrom if chrom.startswith('chr') else 'chr' + str(chrom)
            c2 = c1[3:]
            for c in [c1, c2]:
                self.files[c] = {
                    'root': root,
                    'bed': bed,
                    'bim': bim,
                    'fam': fam,
                }

    def pop_file(self, populations=None, outfile=None):
        """Write temp file with a list of individuals in population."""
        populations = populations if populations else self.individuals.keys()
        if isinstance(populations, str):
            populations = [populations]
        populations = list(populations)
        if outfile and _os.path.isfile(outfile):
            outfile = _os.path.abspath(outfile)
            self.written_files[','.join(populations)] = outfile
            return outfile
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

        if not outfile:
            _, outfile = _temp(prefix='-'.join(populations), dir='/tmp')
        with open(outfile, 'w') as fout:
            fout.write('\n'.join(
                [' '.join(x) for x in zip(pop_ids, pop_ids)]
            ))

        self.written_files[','.join(populations)] = outfile
        return outfile

    def bim_file(self, chrom, snps, outfile=None):
        """Filter a bim file to only include SNPs in snps.

        Parameters
        ----------
        chrom : str
            A chromosome name for the file to work on.
        snps : list_of_tuple
            A list of snps from list_to_rsid_and_locs():
                (rsid, chrom, loc)
        outfile : str, optional
            A file to write to. A tempfile will be created if not is provided.

        Returns
        -------
        file_path : str
            Absolute path to the newly created file.
        """
        if not chrom in self.files:
            raise ValueError(
                '{} is not in our bim files:\n{}'.format(
                    chrom,
                    '\n'.join([i[' bim'] for i in self.files.values()])
                )
            )
        snps = [(str(name), int(loc)) for name, loc in snps]
        if not outfile:
            _, outfile = _temp(
                prefix='filtered_bim.', suffix='.bim', dir='/tmp'
            )
        outfile = _os.path.abspath(outfile)
        with _open_zipped(outfile, 'w') as fout:
            for c, name, cm, pos, a1, a2 in read_bim(self.files[chrom]['bim']):
                if (name, pos) in snps:
                    fout.write('\t'.join([c, name, cm, str(pos), a1, a2]) + '\n')
        return outfile

    def bim_snps(self, chrom, bim_file=None):
        """Return and cache all SNPs in a bim file.

        Note: will only cache rsids in bim if no bim_file is given.

        Parameters
        ----------
        chrom : str
        bim_file : str, optional
            Bim file to work on, otherwise the chromosome is used to pick the
            complete file.

        Returns
        -------
        rsids : frozenset
        """
        if not bim_file:
            if chrom in self.bims:
                with open(self.bims[chrom], 'rb') as fin:
                    return _pickle.load(fin)
            try:
                bim_file = self.files[chrom]['bim']
            except KeyError:
                raise KeyError('{} is not in the bim file list'.format(chrom))
            _, pfl = _temp()
        else:
            pfl = None
        rsids = []
        with open(bim_file) as fin:
            for line in fin:
                if not line.strip():
                    continue
                rsids.append(line.strip().split('\t')[1])
        rsids = frozenset(rsids)
        if pfl:
            with open(pfl, 'wb') as fout:
                _pickle.dump(rsids, fout)
            self.bims[chrom] = pfl
        return rsids

    def many_to_many(self, snps, comp_list, chrom, r2=0.6, populations=None,
                     window_size=50000, fbim_file=None, outfile=None,
                     keep_int_files=False, raise_on_error=False,
                     logfile=_sys.stderr):
        """Get many-to-many LD information using plink.

        Will do a single lookup for all SNPs in snps using plink on a filtered
        bim file generated to only contain the SNPs in comp_list.

        Parameters
        ----------
        snps : list_of_tuple
            list of rsIDs to query in the format:
                (rsid, loc) (from pairs_to_lists())
        comp_list : list_of_tuple
            list of rsIDs to compare to in the format:
                (rsid, loc) (from pairs_to_lists())
        chrom : str
            which chromosome to search
        r2 : float, optional
            r-squared level to use for filtering
        populations : list_of_str, optional
            list of populations to include in the analysis
        window_size : int, optional
            Integer number for window size, default 50kb.
        fbim_file : str, optional
            File to write the filtered bims to, if not provided a temp file
            will be used and deleted.
        outfile : str, optional
            A file root for plink to write to, output will have '.ld' appended
            to it. If not provided, temp file used and deleted after use.
        keep_int_files : bool, optional
            Do not delete intermediate SNP files
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
        window_size = get_length(window_size)
        if not snps or not comp_list:
            raise ValueError('snps and comp_list cannot be null')
        snps = set(snps)
        comp_list = set(comp_list)
        bsnps = self.bim_snps(chrom)
        all_snps = [i for i, _ in snps | comp_list if i in bsnps]
        good  = []
        bad   = []
        rsids = {}
        for rsid, loc in snps:
            if rsid in bsnps:
                good.append(rsid)
                rsids[rsid] = int(loc)
            else:
                bad.append(rsid)
        if bad:
            _sys.stderr.write(
                'The following SNPs are not in the bim file and cannot be ' +
                'queried:\n{}\n'.format(bad)
            )
        del(bsnps)
        # Initialize necessary files
        bfile = self.files[chrom]['root']
        pop_file = self.pop_file(populations)
        if outfile:
            del_file = False
        else:
            _, outfile = _temp(prefix='plink', dir='/tmp')
            _os.remove(outfile)
            del_file = True
        _, snp_file = _temp(prefix='plinksnps', dir='/tmp')
        with open(snp_file, 'w') as fout:
            fout.write('\n'.join(sorted([snp for snp, _ in snps])) + '\n')
        _, comp_file = _temp(prefix='plinkcomp', dir='/tmp')
        with open(comp_file, 'w') as fout:
            fout.write('\n'.join(sorted(all_snps)) + '\n')

        # Build the command
        plink_cmnd = (
            '{plink} --bfile {bfile} --r2 in-phase dprime '
            '--ld-snp-list {snp_file} --extract {comp_file} '
            '--ld-window {window} --keep {ind_file} --out {out}'
        ).format(
            plink=self.plink, bfile=bfile, window=window_size,
            snp_file=snp_file, comp_file=comp_file, ind_file=pop_file,
            out=outfile
        )

        # Run it
        stdout, stderr, code = run(plink_cmnd, raise_on_error)
        # Parse the output file
        if code != 0:
            logfile.write(
                '{}: plink command failed'.format(chrom) +
                'Command: {}\nExit Code: {}\nSTDOUT:\n{}\bSTDERR:\n{}\n'
                .format(plink_cmnd, code, stdout, stderr)
            )
            return {}

        results = {}
        with open(outfile + '.ld') as fin:
            # Check header
            line = fin.readline().strip()
            assert _re.split(r' +', line) == [
                'CHR_A', 'BP_A', 'SNP_A', 'CHR_B', 'BP_B',
                'SNP_B', 'PHASE', 'R2', 'DP'
            ]
            for line in fin:
                f = _re.split(r' +', line.strip())
                snp1, snp2, phased = f[2], f[5], f[6]
                loc1, loc2 = int(f[1]), int(f[4])
                rsquared, dprime = float(f[7]), float(f[8])
                if snp1 not in good:
                    continue
                if snp1 == snp2:
                    continue
                if loc1 != rsids[snp1]:
                    continue
                if rsquared < r2:
                    continue
                if snp1 not in results:
                    results[snp1] = {
                        'chrom': chrom, 'loc': loc1, 'matches': {}
                    }
                try:
                    p1, p2 = phased.split('/')
                    s1a, s2a = p1
                    s1b, s2b = p2
                    lookup = {snp1:  {s1a: s2a, s1b: s2b},
                              snp2: {s2a: s1a, s2b: s1b}}
                except ValueError:
                    lookup = {}
                results[snp1]['matches'][snp2] = {
                    'r2': rsquared, 'dprime': dprime, 'loc': loc2,
                    'phased': phased, 'lookup': lookup
                }

        if del_file:
            pass
            for fl in [_os.path.join('/tmp', f) for f in _os.listdir('/tmp')]:
                if fl.startswith(outfile):
                    _os.remove(fl)
        if not keep_int_files:
            _os.remove(snp_file)
            _os.remove(comp_file)

        return results


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

    If more than 5000 items in rsid list, creates a temp table and joins it.

    Requires about 25 minutes for 250,000 rsIDs.

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
        if len(list1_rs) > 5000:
            conn = _sqlite3.connect(DB_PATH + '/dbsnp{}.db'.format(DB_VERS))
            conn.execute("CREATE TEMP TABLE 'temp' (name)")
            conn.executemany("INSERT INTO temp (name) VALUES (?)", [(i,) for i in list1_rs])
            q = conn.execute("SELECT dbSNP.name, dbSNP.chrom, dbSNP.start, dbSNP.end FROM dbSNP INNER JOIN temp ON dbSNP.name=temp.name").fetchall()
        else:
            db = _dbSNP.DB(DB_PATH, DB_VERS)
            q = []
            for chunk in [list1_rs[i:i + 990] for i in range(0, len(list1_rs), 990)]:
                q += db.query(
                    db.Row.name, db.Row.chrom, db.Row.start, db.Row.end
                ).filter(
                    db.Row.name.in_(chunk)
                ).all()
        rs_lookup = {
            i[0]: (i[0], i[1], i[2]+1) for i in q
            if not i[3]-i[2] > 1
        }

    if list1_locs:
        _sys.stderr.write('Doing position DB lookup.\n')
        db = _dbSNP.DB(DB_PATH, DB_VERS)
        query = {}
        for chrom, loc in [i.split(':') for i in list1_locs]:
            if chrom not in query:
                query[chrom] = []
            query[chrom].append(int(loc)-1)
        loc_lookup = {
            '{}:{}'.format(i.chrom, i.start+1): (i.name, i.chrom, i.start+1)
            for i in db.lookup_locations(query)
            if not i.end-i.start > 1
        }

    failed  = []
    results = []
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
        if s_rs in done:
            continue
        done.append(s_rs)
        results.append((s_rs, s_chrom, s_loc))
    if failed:
        err = '{} not in dbSNP'.format(snp)
        if raise_on_missing:
            raise MissingSNPError(err)
        else:
            _sys.stderr.write(err + '\n')
    return results


def snp_list_to_dict(list1):
    """Convert a list from list_to_rsid_and_locs() to a chrom dict."""
    results = {}
    for s_rs, s_chrom, s_loc in list1:
        if s_chrom not in results:
            results[s_chrom] = []
        results[s_chrom].append((s_rs, int(s_loc)))
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
    if outfile:
        _sys.stderr.write('Writing output.\n')
        if pickle:
            with _open_zipped(outfile, 'wb') as fout:
                _pickle.dump(parsed_list, fout)
        else:
            with _open_zipped(outfile, 'w') as fout:
                for row in parsed_list:
                    rsid, chrom, loc = row
                    fout.write('{}\t{}\t{}\n'.format(rsid, chrom, loc))
    return parsed_list



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
    distance = get_length(distance)
    results = {}
    for chrom, snps in ref_snps.items():
        if chrom not in comp_snps:
            continue
        results[chrom] = distance_filter(snps, comp_snps[chrom], distance)
    return results


def filter_by_ld(pairs, r2=0.6, window_size=50000, populations=None,
                 plink='plink', data_dir=None):
    """Use plink to lookup LD and filter lists by the r2.

    One job per chromosome.

    Parameters
    ----------
    pairs : dict
        The output of filter_by_distance()
    r2 : float
        An r-squared cutoff to use for filtering
    window_size : int
        A distance away from a SNP in list1 to consider SNPs in list2, must be
        a integer of bases
    plink : str, optional
        Path to plink executable, otherwise searches PATH
    data_dir : str, optional
        Path to directory containing 1000genomes data, default set in DATA_DIR
        hardcoded into this file.

    Returns
    -------
    pairs : dict
        snp: {chrom: chrom, loc: loc, matches: {
            snp2: {rsquared: r2, dprime: d', loc: loc2}
    """
    window_size = get_length(window_size)
    plink = PLINK(data=data_dir, plink=plink)
    # Initialize the population file
    plink.pop_file(populations, outfile='.'.join(populations) + '.pops.txt')
    l = 0
    for v in pairs.values():
        for i in [i[2] for i in v if i[2]]:
            l += len(i)
    multi = True if l > 200 else False

    # Run plink jobs
    results = {}
    for chrom in sorted(pairs, key=_chrom_sort):
        data = pairs[chrom]
        snps, comp_list = pairs_to_lists(data)
        results.update(
            plink.many_to_many(
                snps, comp_list, chrom, r2=r2, populations=populations,
                window_size=window_size
            )
        )

    return results


def filter_by_ld_one_by_one(pairs, r2=0.6, populations=None, plink='plink',
                            cluster_jobs=100, force_local=False,
                            **cluster_opts):
    """Use plink to lookup LD and filter lists by the r2 with one job per snp.

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
    cluster_jobs : int, optional
        Number of jobs to run on the cluster if cluster running is enabled, if
        fyrd is not enabled, this ends up being the number of local cores (in
        that case, the max is the result of multiprocessing.cpu_count()).
    force_local : bool, optional
        If True, disallow running jobs on the cluster with fyrd
    cluster_opts : dict, optional
        Optional keyword arguments to pass to fyrd

    """
    plink = PLINK(plink=plink)
    plink.pop_file(populations, outfile=','.join(populations) + '.pops.txt')
    l = 0
    for v in pairs.values():
        for i in [i[2] for i in v if i[2]]:
            l += len(i)
    multi = True if l > 200 else False  # 200 jobs take 50 seconds
    multi = False
    if multi and cluster_jobs > 1:
        # Set up cluster running
        if _fyrd and _fyrd.queue.MODE != 'local' and not force_local:
            print('Running jobs on the cluster with fyrd.')
            fyrd = True
        else:
            print('Running jobs in parallel locally')
            cores = max(cluster_jobs, _mp.cpu_count())
            pool = _mp.Pool(cores)
            fyrd = False
        # Bundle jobs
        bundle_count = int(l/cluster_jobs)
        bundle_count = max(bundle_count, MIN_BUNDLE)
        count = bundle_count
        # Get walltime based on 6s per job + 40 minutes
        # Average on a standard machine is 2 seconds, so this is
        # an overestimate
        time = _to_time((bundle_count*6))
        _sys.stderr.write('Estimated max time per job {} '.format(time) +
                          '(walltime req will add 40 min to this)\n')

    jobs = [[]]
    n = 0
    skipped = 0
    for chrom in sorted(pairs, key=_chrom_sort):
        snps = pairs[chrom]
        args = (chrom, r2, populations, False)
        # Bundle jobs
        if not multi:
            # Bundle by chromosome if not multithreading
            bundle_count = len(snps)
            count = bundle_count
        for snp, _, comp_list in snps:
            if not comp_list:
                skipped += 1
                continue
            comp_snps = [rsid for rsid, loc in comp_list]
            if not count:
                count = bundle_count
                jobs.append([])
                n += 1
            jobs[n].append((snp, comp_snps) + args)
            count -= 1

    if skipped:
        _sys.stderr.write('{} jobs had no comparison list and were skipped\n'
                          .format(skipped))
    # Submit jobs
    running = []
    results = {}
    for jobset in jobs:
        if multi:
            if fyrd:
                time = _to_time((len(jobset)*6)+2400)
                running.append(
                    _fyrd.submit(
                        _ld_job, (jobset, plink),
                        walltime=time, mem='8G', cores=4, **cluster_opts)
                )
            else:
                running.append(pool.apply_async(_ld_job, jobset))
        else:
            # Just run jobs per chromosome
            results.update(_ld_job(jobset, plink))

    if multi:
        print('Getting results from parallel jobs')
        if fyrd:
            todo = len(running)
            while todo:
                for job in running:
                    if job.done:
                        results.update(job.get())
                        todo -= 1
                _sleep(2)
        else:
            for job in running:
                results.update(job.get())
            pool.close()

    # Add all data back again
    filtered = {}
    for chrom, data in pairs.items():
        for snp, loc, comp_list in data:
            if snp not in results:
                continue
            if results[snp]:
                filtered[snp] = {'loc': loc, 'chrom': chrom}
                filtered[snp]['matches'] = {}
                for snp2, pos2 in comp_list:
                    if snp2 in results[snp]:
                        filtered[snp]['matches'][snp2] = results[snp][snp2]
                        filtered[snp]['matches'][snp2]['loc'] = pos2
    return filtered


def _ld_job(jobs, plink):
    """Private function to run plink, see filter_by_ld().

    Parameters
    ----------
    jobs : list_of_tuple
        A list of job arguments:
            (snp, comp_snps, chrom, r2, populations, raise_on_error, logfile)
    plink : PLINK
        An initialized plink class
    """
    snp_results = {}
    for job in jobs:
        snp_results[job[0]] = plink.one_to_many(*job)
    return snp_results


def _to_time(s):
    """Convert seconds to 00:00:00 format."""
    m, s = divmod(s, 60)
    h, m = divmod(m, 60)
    return '{}:{:02}:{:02}'.format(h, m, s)


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
                   plink='plink', lists_preparsed=False, one_to_many=False,
                   dbsnp_loc=None, dbsnp_ver=None, data_dir=None,
                   return_dataframe=False, return_snp_pairs=False,
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
    lists_preparsed : bool, optional
        If the input lists are already in three-tuples will all data acounted
        for, skip the dbSNP phase.
    one_to_many : bool, optional
        Force plink to run one job per lookup SNP, much slower.
    dbsnp_loc : str, optional
        Path to dbSNP sqlite database files
    dbsnp_ver : int : optional
        Version of dbSNP to use
    data_dir : str, optional
        Path to 1000genomes files
    return_dataframe : bool, optional
        Return a dataframe instead of JSON
    return_snp_pairs : bool, optional
        Return a dictionary of SNP_Link objects
    raise_on_error : bool, optional
        Raise an Exception if there is a problem with any of the SNPs in the
        list, otherwise errors are silently ignored.

    Returns
    -------
    JSON : dict
        If return_dataframe and return_snp_pairss are False, a json-style
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
    distance = get_length(distance)
    if not lists_preparsed:
        print('Getting SNP info for list1')
        list1 = list_to_rsid_and_locs(list1, raise_on_error)
        print('Getting SNP info for list2')
        list2 = list_to_rsid_and_locs(list2, raise_on_error)

    snps1 = snp_list_to_dict(list1)
    snps2 = snp_list_to_dict(list2)

    print('Getting all possible pairs within {}'.format(distance))
    pairs = filter_by_distance(snps1, snps2, distance=distance)
    l = 0
    for v in pairs.values():
        for i in [i[2] for i in v]:
            l += len(i)
    print('Running plink filters on {} calculations'.format(l))
    if one_to_many:
        # Not a good method, but kept for comparisons
        pairs = filter_by_ld_one_by_one(
            pairs, r2, populations=populations, plink=plink,
            force_local=True
        )
    else:
        pairs = filter_by_ld(pairs, r2, populations=populations, plink=plink,
                             window_size=distance, data_dir=data_dir)
    print('Done')
    if return_dataframe:
        return pairs_to_dataframe(pairs)
    if return_snp_pairs:
        return pairs_to_SNP_Pairs(pairs, populations=populations)
    return pairs


###############################################################################
#                              Helper Functions                               #
###############################################################################


class MissingSNPError(Exception):

    """Exceception to catch missing SNPs."""

    pass


class BadSNPError(Exception):

    """Exceception to catch missing SNPs."""

    pass


def pairs_to_lists(pairs):
    """Convert list of pairs into to two tuple lists.

    Removes any SNPs with an empty match list.

    Parameters
    ----------
    pairs : list
        From filter_by_distance():
            (rsid, loc, set((rsid, loc)))

    Returns
    -------
    snp_list : list_of_tuple
    comp_list : list_of_tuple
        (rsid, loc)
    """
    query = []
    comp  = []
    for snp, loc, matches in pairs:
        if not matches:
            continue
        query.append((snp, loc))
        for snp2, loc2 in matches:
            comp.append((snp2, loc2))
    return set(query), set(comp)


def get_length(length_string):
    """Convert a string length like 50kb to an integer length."""
    # Calculate length_string
    if isinstance(length_string, str):
        if length_string.isdigit():
            return int(length_string)
        dist, mod = _re.match(r'(\d+)(\D+)', length_string).groups()
        dist = int(dist)
        mod = mod.lower()
        if mod:
            if not mod in ['kb', 'mb']:
                raise ValueError('Cannot parse {}, must be in kb or mb only'
                                 .format(length_string))
            if mod == 'kb':
                dist = dist * 1000
            elif mod == 'mb':
                dist = dist * 1000000
    elif isinstance(length_string, int):
        return length_string
    else:
        raise ValueError('length_string must be either a string or an integer, is '
                         '{}'.format(type(length_string)))
    return dist


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


def read_bim(bim_file):
    """Yields a tuple for each line in bim_file.

    Parameters
    ----------
    bim_file : str
        Path to a bim file, can be zipped or an open file handle.

    Yields
    ------
    chromsome : str
    name : str
    cm : str
        Position in centimorgans, usually 0
    position : int
    allele_1 : str
    allele_2 : str
    """
    with _open_zipped(bim_file) as fin:
        for line in fin:
            chrom, name, cm, loc, a1, a2 = line.split('\t')
            yield chrom, name, cm, int(loc), a1, a2


def _open_zipped(infile, mode='r'):
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


def argument_parser():
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

    # Optional files
    paths = analyze.add_argument_group(
        'data_paths',
        'Paths to 1000genomes and dbSNP'
    )
    paths.add_argument('--dbsnp-path', default=DB_PATH,
                       help='Path to dbSNP directory')
    paths.add_argument('--dbsnp-ver', default=DB_VERS,
                       help='dbSNP version to use')
    paths.add_argument('--1000genomes', default=DATA_DIR,
                       help='1000genomes data files')

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

    parser = argument_parser()

    args = parser.parse_args(argv)

    if not args.mode:
        parser.print_usage()
        _sys.stdout.write(general_usage)
        _sys.stderr.write('Mode required.\n')
        return 2

    with _open_zipped(args.snp_list) as fin:
        snp_list = fin.read().strip().split('\n')

    if args.mode == 'preprocess':
        outfile = args.outfile if args.outfile else _sys.stdout
        save_list(snp_list, outfile)
        return 0

    if not args.outfile and args.pandas:
        _sys.stderr.write('Cannot write pandas DataFrame to STDOUT\n')
        return 1

    with _open_zipped(args.comp_list) as fin:
        comp_list = fin.read().strip().split('\n')

    if args.populations:
        pops = args.populations.split(',')
    else:
        pops = None

    out = comp_snp_lists(
        snp_list, comp_list, pops,
        r2=float(args.r2),
        distance=args.distance,
        return_dataframe=True,
    )

    if args.pandas:
        out.to_pickle(args.outfile)
    else:
        out.to_csv(args.outfile, sep='\t')

    return 0

if __name__ == '__main__' and '__file__' in globals():
    _sys.exit(main())
