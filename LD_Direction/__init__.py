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
Version: 0.1a

LDLink
------
https://analysistools.nci.nih.gov/LDlink
https://github.com/CBIIT/nci-webtools-dceg-linkage

This tool uses the above API entirely, thanks to Mitchell Machiela and the
Chanock lab for writing that tool!

Citation
~~~~~~~~
Machiela MJ, Chanock SJ. LDlink a web-based application for exploring
population-specific haplotype structure and linking correlated alleles of
possible functional variants. Bioinformatics. 2015 Jul 2. PMID: 26139635.
"""
__version__ = '0.1a'
