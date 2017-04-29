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
import sys as _sys
import argparse as _argparse


def comp_snp_lists(list1, list2, populations=None, r2=0.9, distance='50kb',
                   output_table=False, return_dataframe=False):
    """Compare two SNP lists."""
    pass


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
