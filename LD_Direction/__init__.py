# -*- coding: utf-8 -*-
"""
A simple package to find the Linkage Disequilibrium between any two SNP lists.

The ultimate idea is to be able to pull out all SNPs in list two that are
within some distance of a SNP in list1 (e.g. 50kb in either direction), and
ask what the LD between them is. For all those that meet some R-squared cutoff,
we will return some LD-stats (R-squared, D-prime, p-value) and the direction
of agreement, along with counts for agreements.

Info
----
Author: Michael D Dacre, mike.dacre@gmail.com
Original author (LDLink): Mitchell Machiela
Organization: Stanford University
License: MIT License, property of Stanford and NCI, use as you wish
Created: 2017-21-21 10:04
Version: 0.2.0b1

Modules
-------
compare_snp_lists
    provides functions to do a one-to-many comparison for arbitrarily large
    lists of SNPs (limited only by memory)
snp_link
    provides functions and classes to do one-to-one LD comparisons on either
    just two SNPs, or a list of SNP pairs
LDpair
    Adapted from LDlink, provides LD calculation tools directly from the VCF
    file without the need for plink

PLINK
-----
https://www.cog-genomics.org/plink
https://github.com/chrchang/plink-ng/

compare_snp_lists depends entirely upon plink and could not function without
it. Many, many thanks to Christopher Chang for creating it, and for his help
with getting this package working.

LDLink
------
https://analysistools.nci.nih.gov/LDlink
https://github.com/CBIIT/nci-webtools-dceg-linkage

This package makes extensive use of the LDlink API entirely, thanks to Mitchell
Machiela and the Chanock lab for writing that tool!

Citations
---------
Chang, Christopher C and Chow, Carson C and Tellier, Laurent C A M and
    Vattikuti, Shashaank and Purcell, Shaun M and Lee, James J. Second-
    generation {PLINK}: {R}ising to the challenge of larger and richer
    datasets. GigaScience. 2015. DOI 10.1186/s13742-015-0047-8
Machiela MJ, Chanock SJ. LDlink a web-based application for exploring
    population-specific haplotype structure and linking correlated alleles of
    possible functional variants. Bioinformatics. 2015 Jul 2. PMID: 26139635.
"""
__version__ = '0.2.0b2'

__all__ = ['snp_link', 'compare_snp_lists']

# Make core functions easily available
from . import compare_snp_lists
from . import snp_link
