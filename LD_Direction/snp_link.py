#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Investigate linkage disequilibrium between pairs of SNPs.

Uses the LDLink API to search for pairs of SNPs and builds a simple class
(SNP_Pair) with the results.

Examples
--------

SNPs in LD
~~~~~~~~~~
>>>snp = compare_two_variants('rs11121194', 'rs1884352', ['ESN'])

>>>snp
SNP_Pair<rs11121194(['C', 'T']:rs1884352(['G', 'A']) R2: 0.9222 (P: <0.0001)>

>>>snp.rsquared
0.922

>>>snp.lookup_other(1, 'C')
'G'

>>>snp.lookup_other('rs11121194', 'C')
'G'

>>>snp.inld
True

SNPs not in LD
~~~~~~~~~~~~~~
>>>snp2 = compare_two_variants('rs75830605', 'rs77977823', ['CEU'])

>>>snp2.inld
False

>>>snp2.lookup
{}

>>>snp2.lookup_other('rs77977823', 'C') # Prints error message to STDERR
None
===============================================================================
"""
import re as _re
import sys as _sys
from time import sleep as _sleep
from urllib.request import urlopen as _get

import pandas as _pd
from numpy import nan as _nan

from . import LDpair as _ld

SLEEP_TIME = 1.0

###############################################################################
#                          Class to hold the linkage                          #
###############################################################################


class SNP_Pair(object):

    """Association information from a pair of SNPs, created from LDpair.

    Attributes
    ----------
    inld : bool
        True if SNPs are in LD, False if SNPs are in Linkage Equilibrium.
    snp1 : str
    snp2 : str
    chrom :  str
    loc1 : int
        Position of snp1 on chromosome
    loc2 : int
        Position of snp2 on chromosome
    dprime : float_or_nan
    rsquared : float_or_nan
    chisq : float_or_nan
    p : float_or_nan
    p_str : str
        String representation of the p-value
    populations : list
        List of 1000genomes populations
    table : pandas.DataFrame
        Pandas DataFrame of allele counts for both SNPs in the above
        populations. Rows are snp1, columns are snp2. Multiindexed.
    alleles : dict
        Dictionary of alleles by snp
    lookup : dict
        Dictionary of allele in other SNP given an allele in one SNP.

    Methods
    -------
    lookup_other(self, snp, allele)
        Get the allele of a SNP given an allele in the other SNP.

    """

    def __init__(self, x):
        """Parse an input dictionary from LDpair."""
        self.input_dict = x

        self.snp1 = x['snp1']['name']
        self.snp2 = x['snp2']['name']

        self.chrom, self.loc1 = x['snp1']['coord'].split(':')
        c2, self.loc2 = x['snp2']['coord'].split(':')
        if c2 != self.chrom:
            _sys.stderr.write('{}: WARNING: SNPs on different chromosomes.'
                              .format([self.snp1, self.snp2]))
            self.chrom = 'multiple'
        self.loc1 = int(self.loc1)
        self.loc2 = int(self.loc2)

        self.snp1_rs = x['snp1']['rsnum']
        self.snp2_rs = x['snp2']['rsnum']

        self.haplotypes = x['haplotypes']

        self.populations = x['populations']

        if x['corr_alleles'][0].endswith('are in linkage equilibrium'):
            self.lookup = {}
            self.inld = False
        else:
            self.lookup = x['corr_dictionary']
            self.inld = True

        self.alleles = {
            self.snp1: [
                x['snp1']['allele_1']['allele'],
                x['snp1']['allele_2']['allele']
            ],
            self.snp2 : [
                x['snp2']['allele_1']['allele'],
                x['snp2']['allele_2']['allele']
            ]
        }

        # Make a pandas table of haplotypes
        cols = _pd.MultiIndex.from_tuples(
            [(self.snp2, i) for i in self.alleles[self.snp1]]
        )
        indx = _pd.MultiIndex.from_tuples(
            [(self.snp1, i) for i in self.alleles[self.snp2]]
        )
        t = x['two_by_two']['cells']
        self.table = _pd.DataFrame(
            [[int(t['c11']), int(t['c12'])],
             [int(t['c21']), int(t['c22'])]],
            columns=cols, index=indx
        )

        self.dprime = x['statistics']['d_prime']
        try:
            self.dprime = float(self.dprime)
        except ValueError:
            self.dprime = _nan

        self.rsquared = x['statistics']['r2']
        try:
            self.rsquared = float(self.dprime)
        except ValueError:
            self.rsquared = _nan

        self.chisq = x['statistics']['chisq']
        try:
            self.chisq = float(self.dprime)
        except ValueError:
            self.chisq = _nan

        p = x['statistics']['p']
        try:
            self.p = 0.0 if p == '<0.0001' else float(p)
        except ValueError:
            self.p = _nan
        self.p_str = p

    def lookup_other(self, snp, allele):
        """Return the linked allele for a given snp.

        Parameters
        ----------
        snp : int_or_str
            Either 1, 2 for SNP 1/2 or rsID.
        allele :str
            The allele for snp.

        Returns
        -------
        str_or_None
            Linked allele for other SNP, unless not in LD.

        Raises
        ------
        ValueError
            If snp or allele contain non-permitted values.
        """
        if not self.inld:
            print('Cannot lookup allele, SNPs in linkage equilibrium.',
                  file=_sys.stderr)
            return

        if isinstance(snp, int):
            if snp == 1:
                snp = self.snp1
            elif snp == 2:
                snp = self.snp2
            else:
                raise ValueError('Invalid value for SNP')
        if snp not in [self.snp1, self.snp2]:
            raise ValueError('SNP must be one of {}'.format([self.snp1, self.snp2]))

        allele = allele.upper()

        if allele not in self.alleles[snp]:
            raise ValueError('Allele {} is invalid for SNP {}, possible values are {}'
                             .format(allele, snp, self.alleles[snp]))

        return self.lookup[snp][allele]

    def __repr__(self):
        """Return infomation"""
        return "SNP_Pair<{snp1}({snp1_a}):{snp2}({snp2_a}) R2: {r2} (P: {P})>".format(
            snp1=self.snp1, snp1_a=self.alleles[self.snp1],
            snp2=self.snp2, snp2_a=self.alleles[self.snp2],
            r2=self.rsquared, P=self.p_str
        )

    def __str__(self):
        """Print summary"""
        return _ld.print_summary(self.input_dict)


###############################################################################
#                                Main Functions                               #
###############################################################################


def compare_variants(snp_list: list, populations: list, rsquared: float=None):
    """Yield SNP_Pair objects for every pair of variants in snp_list.

    Will wait some amount of time between calls, defined by SLEEP_TIME.

    Parameters
    ----------
    snp_list : list of tuples
        Format: [(var1, var2), ...]

    populations : list
        List of 1000genomes populations (e.g. YRI, ESN), string is converted to
        a list (allowing single population lookup).

    rsquared : float
        An R-squared cutoff to use for returning results, those with an
        r-squared below that value will not be returned. Note that as the
        query must be executed prior to checking the r-squared, a large number
        of requests with low r-squared values will result in long pauses prior
        to a value being yielded.

    Yields
    ------
    SNP_Pair
        list of SNP_Pair objects for each pair of variants in snp_list.
    """
    assert isinstance(snp_list, list)
    assert isinstance(snp_list[0], tuple)
    assert len(snp_list[0]) == 2
    for snp1, snp2 in snp_list:
        snp_pair = compare_two_variants(snp1, snp2, populations)
        if rsquared and isinstance(rsquared, float):
            if snp_pair.rsquared < rsquared:
                print('{} failed r-squared cutoff'.format(repr(snp_pair)))
                _sleep(SLEEP_TIME)
                continue
        yield snp_pair
        _sleep(SLEEP_TIME)


def compare_two_variants(var1: str, var2: str, populations: list) -> SNP_Pair:
    """Return a SNP_pair class for any two rsids.

    Uses the LDpair API:
        https://analysistools.nci.nih.gov/LDlink/?tab=ldpair

    Parameters
    ----------
    var1/var2 : str
        rsIDs to compare
    populations : list
        list of 1000genomes populations (e.g. YRI, ESN), string is converted to
        a list (allowing single population lookup).

    Returns
    -------
    SNP_Pair
        A SNP_Pair class with information about SNP linkage.
    """
    assert isinstance(var1, str)
    assert isinstance(var2, str)
    if isinstance(populations, str):
        populations = [populations]

    return SNP_Pair(_ld.calculate_pair(var1, var2, populations))


###############################################################################
#                                   Helpers                                   #
###############################################################################


# Parse allele linkage data from last two lines.
correlation_lookup = _re.compile(
    r'^(rs[0-9]+)\(([ATGC])\) allele is correlated with ' +
    r'(rs[0-9]+)\(([ATGC])\) allele$'
)


def get_snps(x: str) -> tuple:
    """Parse a SNP line and return name, chromsome, position."""
    snp, loc = x.split(' ')
    chrom, position = loc.strip('()').split(':')
    return snp, chrom, int(position)
