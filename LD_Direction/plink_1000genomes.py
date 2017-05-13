#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Methods to run plink commands on the 1000genomes phase 3 dataset.
"""
import os as _os
import re as _re
import sys as _sys
import pickle as _pickle
from tempfile import mkstemp as _temp

from . import _run

# Set data directories
DATA_DIR = "/godot/1000genomes/1000GP_Phase3"

POPULATIONS = ["ALL", "AFR", "AMR", "EAS", "EUR", "SAS", "ACB", "ASW", "BEB",
               "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH",
               "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL",
               "PJL", "PUR", "STU", "TSI", "YRI"]

__all__ = ["PLINK", "read_bim"]

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
        with _run.open_zipped(outfile, 'w') as fout:
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
                     window_size=50000, outfile=None,
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
        window_size = _run.get_length(window_size)
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
        stdout, stderr, code = _run.run(plink_cmnd, raise_on_error)
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
        stdout, stderr, code = _run.run(plink_cmnd, raise_on_error)
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
#                              Helper Functions                               #
###############################################################################


class MissingSNPError(Exception):

    """Exceception to catch missing SNPs."""

    pass


class BadSNPError(Exception):

    """Exceception to catch missing SNPs."""

    pass


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
    with _run.open_zipped(bim_file) as fin:
        for line in fin:
            chrom, name, cm, loc, a1, a2 = line.split('\t')
            yield chrom, name, cm, int(loc), a1, a2
