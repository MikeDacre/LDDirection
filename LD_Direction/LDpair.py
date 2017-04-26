#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Modified LDpair from:

https://github.com/CBIIT/nci-webtools-dceg-linkage/blob/master/LDlink/LDpair.py
"""
import os
import sys
import json
import math
import time
import sqlite3
import argparse
import subprocess

DB_SNP_VER = '149'

POPULATIONS = ["ALL", "AFR", "AMR", "EAS", "EUR", "SAS", "ACB", "ASW", "BEB",
               "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH",
               "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL",
               "PJL", "PUR", "STU", "TSI", "YRI"]

# Set data directories
DATA_DIR = "/godot/1000genomes/1000GP_Phase3"
SNP_DB   = "/godot/dbsnp/dbsnp144.db"


class SNP_Lookup_Failure(Exception):

    """Simple Exception."""

    pass


def output_json(output):
    """Return JSON formated string."""
    return json.dumps(output, sort_keys=True, indent=2)


def run_cmnd(cmnd):
    """Run a command and return the output split into a list by newline."""
    return subprocess.check_output(
        cmnd, shell=True
    ).decode().strip().split('\n')


def get_snp_info(snp):
    """Use dbSNP and 1000genomes tabix to get SNP information.

    Returns
    -------
    snp_chrom : str
    snp_loc : int
    geno : list
        Line from a VCF file corresponding to this SNP split on tab
    allele : dict
        Dictionary with numbered genotype to allele data, e.g.:
            {'.': ['.', '.'], './.': ['.', '.'],
             '0': ['G', '.'], '0|0': ['G', 'G'], '0|1': ['G', 'A'],
             '1': ['A', '.'], '1|0': ['A', 'G'], '1|1': ['A', 'A']}
    head : list
        The headers from the 1000genomes VCF file
    allele1 : str
    allele2 : str

    Raises
    ------
    SNP_Lookup_Failure
        Any failure at any step (getting location from rsID, lookup in VCF file,
        etc) will raise a SNP_Lookup_Failure. Any other Exception should not
        happen.
    """
    # Find RS numbers in dbSNP
    if snp.startswith('rs'):
        # Connect to snp database
        conn = sqlite3.connect(SNP_DB)
        conn.text_factory = str
        cur = conn.cursor()

        id = snp.strip("rs")
        t = (id,)
        cur.execute("SELECT * FROM tbl_" + id[-1] + " WHERE id=?", t)
        snp_coord = cur.fetchone()
        if not snp_coord:
            raise SNP_Lookup_Failure('{} is not in dbSNP build {}.'
                                     .format(snp, DB_SNP_VER))
        snp_chrom = d[1]
        snp_loc = d[2]
    elif ':' in snp:
        snp_chrom, snp_loc = snp.split(':')
    else:
        raise SNP_Lookup_Failure(snp + ' is malformed, does not appear to be '
                                 'rsID or location (chr:position)')

    snp_chrom = snp_chrom if snp_chrom.startswith('chr') else 'chr' + str(snp_chrom)

    # Get SNP info from 1000genomes
    vcf_file = (
        DATA_DIR +
        "/ALL.{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
        .format(snp_chrom)
    )
    # Strip 'chr' from chromosome
    tabix_cmd = "tabix {0} {1}:{2}-{2} | grep -v -e END".format(
        vcf_file, snp_chrom[3:], snp_loc
    )
    vcf = run_cmnd(tabix_cmd)

    # Import SNP VCF files
    if len(vcf) == 0:
        raise SNP_Lookup_Failure(snp + " is not in 1000G reference panel.")
    elif len(vcf) == 1:
        geno = vcf[0].strip().split()
    else:
        geno = []
        for line in vcf:
            l = line.strip().split()
            if l[2] == snp:
                geno = l
                break
        if geno == []:
            raise SNP_Lookup_Failure(snp + " is not in 1000G reference panel.")

    if "," in geno[3] or "," in geno[4]:
        raise SNP_Lookup_Failure(snp + " is not a biallelic variant.")

    if len(geno[3]) == 1 and len(geno[4]) == 1:
        snp_a1 = geno[3]
        snp_a2 = geno[4]
    elif len(geno[3]) == 1 and len(geno[4]) > 1:
        snp_a1 = "-"
        snp_a2 = geno[4][1:]
    elif len(geno[3]) > 1 and len(geno[4]) == 1:
        snp_a1 = geno[3][1:]
        snp_a2 = "-"
    elif len(geno[3]) > 1 and len(geno[4]) > 1:
        snp_a1 = geno[3][1:]
        snp_a2 = geno[4][1:]

    allele = {
        "0|0": [snp_a1, snp_a1],
        "0|1": [snp_a1, snp_a2],
        "1|0": [snp_a2, snp_a1],
        "1|1": [snp_a2, snp_a2],
        "0":   [snp_a1, "."],
        "1":   [snp_a2, "."],
        "./.": [".", "."],
        ".":   [".", "."]
    }

    # Get headers
    tabix_header_cmd = "tabix -H {} | grep CHROM".format(vcf_file)
    head = run_cmnd(tabix_header_cmd)[0].strip().split()

    return snp_chrom, snp_loc, geno, allele, head, snp_a1, snp_a2


def calculate_pair(snp1, snp2, pops, write_summary=None, return_json=False):
    """Find LD information for any two SNPs.

    Parameters
    ----------
    snp1/snp2 : str
        rsIDs to compare
    pops : list
        list of populations (e.g. YRI, ESN).
    write_summary : str
        a file to write output to. if 'STDOUT' or 'STDERR' will print to those
        locations
    return_json : bool
        if True, return formatted JSON instead of a dictionary, this has
        the side-effect of returning a JSON error string instead of raising an
        exception on error.

    Returns
    -------
    dict_or_str
        A dictionary of structured data or formatted JSON.
    """
    # Create JSON output
    output = {}

    # Get SNP info
    try:
        (snp1_chrom, snp1_loc, geno1,
         allele1, head1, snp1_a1, snp1_a2) = get_snp_info(snp1)
        (snp2_chrom, snp2_loc, geno2,
         allele2, head2, snp2_a1, snp2_a2) = get_snp_info(snp2)
    except SNP_Lookup_Failure as e:
        if return_json:
            output["error"] = str(e)
            return output_json(output)
        else:
            raise

    # Check if SNPs are on the same chromosome
    if snp1_chrom != snp2_chrom:
        output["warning"] = ("{} and {} are on different chromosomes"
                             .format(snp1, snp2))

    for geno, snp in [(geno1, snp1), (geno2, snp2)]:
        if not snp.startswith('rs'):
            continue
        if geno[2] != snp:
            if "warning" not in output:
                output["warning"] = ''
            else:
                output["warning"] += '. '
            output["warning"] += (
                "Genomic position for query variant1 ({}) does not match RS "
                "number at 1000G position ({})"
            ).format(snp, geno[2])
    snp1_rs = snp1 if snp1.startswith('rs') else geno1[2]
    snp2_rs = snp2 if snp2.startswith('rs') else geno2[2]

    if geno1[1] != str(snp1_loc):
        err = "VCF File does not match variant coordinates for SNP1."
        if return_json:
            output["error"] = err
            return output_json(output)
        else:
            raise SNP_Lookup_Failure(err)
    if geno2[1] != str(snp2_loc):
        err = "VCF File does not match variant coordinates for SNP2."
        if return_json:
            output["error"] = err
            return output_json(output)
        else:
            raise SNP_Lookup_Failure(err)

    ### Combine phased genotypes ###
    geno = {}
    for i in range(9, len(head1)):
        geno[head1[i]] = [allele1[geno1[i]], ".."]

    for i in range(9, len(head2)):
        if head2[i] in geno:
            geno[head2[i]][1] = allele2[geno2[i]]

    # Select desired ancestral populations
    pop_file = os.path.join(DATA_DIR,
                            'integrated_call_samples_v3.20130502.ALL.panel')
    individuals = {}
    with open(pop_file) as fin:
        assert fin.readline() == 'sample\tpop\tsuper_pop\tgender\t\t\n'
        for line in fin:
            ind, pop, _, _ = line.split('\t')
            if pop not in individuals:
                individuals[pop] = []
            individuals[pop].append(ind)

    pop_ids  = []
    bad_pops = []
    for pop_i in pops:
        if pop_i in individuals:
            pop_ids += individuals[pop_i]
        else:
            bad_pops.append(pop_i)
    if bad_pops:
        err = (
            "{} are not ancestral populations. Choose one of the following "
            "ancestral populations: {}"
        ).format(bad_pops, POPULATIONS)
        if return_json:
            output["error"] = err
            return output_json(output)
        else:
            raise SNP_Lookup_Failure(err)

    pop_ids = list(set(pop_ids))

    # Extract haplotypes
    hap = {}
    for ind in pop_ids:
        if ind in geno:
            hap1 = geno[ind][0][0] + "_" + geno[ind][1][0]
            hap2 = geno[ind][0][1] + "_" + geno[ind][1][1]

            if hap1 in hap:
                hap[hap1] += 1
            else:
                hap[hap1] = 1

            if hap2 in hap:
                hap[hap2] += 1
            else:
                hap[hap2] = 1

    # Remove missing haplotypes
    keys = list(hap.keys())
    for key in keys:
        if "." in key:
            hap.pop(key, None)

    # Check all haplotypes are present
    if len(hap) != 4:
        snp1_a = [snp1_a1, snp1_a2]
        snp2_a = [snp2_a1, snp2_a2]
        haps = [snp1_a[0] + "_" + snp2_a[0], snp1_a[0] + "_" + snp2_a[1],
                snp1_a[1] + "_" + snp2_a[0], snp1_a[1] + "_" + snp2_a[1]]
        for i in haps:
            if i not in hap:
                hap[i] = 0

    # Sort haplotypes
    A = hap[sorted(hap)[0]]
    B = hap[sorted(hap)[1]]
    C = hap[sorted(hap)[2]]
    D = hap[sorted(hap)[3]]
    N = A + B + C + D
    tmax = max(A, B, C, D)

    hap1 = sorted(hap, key=hap.get, reverse=True)[0]
    hap2 = sorted(hap, key=hap.get, reverse=True)[1]
    hap3 = sorted(hap, key=hap.get, reverse=True)[2]
    hap4 = sorted(hap, key=hap.get, reverse=True)[3]

    delta = float(A * D - B * C)
    Ms = float((A + C) * (B + D) * (A + B) * (C + D))
    if Ms != 0:
        # D prime
        if delta < 0:
            D_prime = abs(delta / min((A + C) * (A + B), (B + D) * (C + D)))
        else:
            D_prime = abs(delta / min((A + C) * (C + D), (A + B) * (B + D)))

        # R2
        r2 = (delta**2) / Ms

        # P-value
        num = (A + B + C + D) * (A * D - B * C)**2
        denom = Ms
        chisq = num / denom
        p = 2 * (1 - (0.5 * (1 + math.erf(chisq**0.5 / 2**0.5))))

    else:
        D_prime = "NA"
        r2 = "NA"
        chisq = "NA"
        p = "NA"

    # Find Correlated Alleles
    if r2 != "NA" and r2 > 0.1:

        # Expected Cell Counts
        eA = (A + B) * (A + C) / N
        eB = (B + A) * (B + D) / N
        eC = (C + A) * (C + D) / N
        eD = (D + C) * (D + B) / N

        # Calculate Deltas
        dA = (A - eA)**2
        dB = (B - eB)**2
        dC = (C - eC)**2
        dD = (D - eD)**2
        dmax = max(dA, dB, dC, dD)

        outstr = "{} ({}) allele is correlated with {} ({})"
        corr_lookup = {}
        if dA == dB == dC == dD:
            if tmax == A or tmax == D:
                corr_lookup[snp1] = {
                    sorted(hap)[0].split("_")[0]: sorted(hap)[0].split("_")[1],
                    sorted(hap)[2].split("_")[0]: sorted(hap)[1].split("_")[1]
                }
                corr_lookup[snp2] = {
                    sorted(hap)[0].split("_")[1]: sorted(hap)[0].split("_")[0],
                    sorted(hap)[1].split("_")[1]: sorted(hap)[2].split("_")[0]
                }
            else:
                corr_lookup[snp1] = {
                    sorted(hap)[0].split("_")[0]: sorted(hap)[1].split("_")[1],
                    sorted(hap)[2].split("_")[0]: sorted(hap)[0].split("_")[1]
                }
                corr_lookup[snp2] = {
                    sorted(hap)[1].split("_")[1]: sorted(hap)[0].split("_")[0],
                    sorted(hap)[0].split("_")[1]: sorted(hap)[2].split("_")[0]
                }
        elif dmax == dA or dmax == dD:
            corr_lookup[snp1] = {
                sorted(hap)[0].split("_")[0]: sorted(hap)[0].split("_")[1],
                sorted(hap)[2].split("_")[0]: sorted(hap)[1].split("_")[1]
            }
            corr_lookup[snp2] = {
                sorted(hap)[0].split("_")[1]: sorted(hap)[0].split("_")[0],
                sorted(hap)[1].split("_")[1]: sorted(hap)[2].split("_")[0]
            }
        else:
            corr_lookup[snp1] = {
                sorted(hap)[0].split("_")[0]: sorted(hap)[1].split("_")[1],
                sorted(hap)[2].split("_")[0]: sorted(hap)[0].split("_")[1]
            }
            corr_lookup[snp2] = {
                sorted(hap)[1].split("_")[1]: sorted(hap)[0].split("_")[0],
                sorted(hap)[0].split("_")[1]: sorted(hap)[2].split("_")[0]
            }
        corr_alleles = [
            outstr.format(snp1, a1, snp2, a2)
            for a1, a2 in corr_lookup[snp1].items()
        ]
    else:
        corr_alleles = [snp1 + " and " + snp2 + " are in linkage equilibrium"]

    # Create dictionary of all data
    snp_1 = {}
    snp_1["rsnum"] = snp1_rs
    snp_1["coord"] = '{}:{}'.format(snp1_chrom, snp1_loc)

    snp_1_allele_1 = {}
    snp_1_allele_1["allele"] = sorted(hap)[0].split("_")[0]
    snp_1_allele_1["count"] = str(A + B)
    snp_1_allele_1["frequency"] = str(round(float(A + B) / N, 3))
    snp_1["allele_1"] = snp_1_allele_1

    snp_1_allele_2 = {}
    snp_1_allele_2["allele"] = sorted(hap)[2].split("_")[0]
    snp_1_allele_2["count"] = str(C + D)
    snp_1_allele_2["frequency"] = str(round(float(C + D) / N, 3))
    snp_1["allele_2"] = snp_1_allele_2
    output["snp1"] = snp_1

    snp_2 = {}
    snp_2["rsnum"] = snp2_rs
    snp_2["coord"] = '{}:{}'.format(snp2_chrom, snp2_loc)

    snp_2_allele_1 = {}
    snp_2_allele_1["allele"] = sorted(hap)[0].split("_")[1]
    snp_2_allele_1["count"] = str(A + C)
    snp_2_allele_1["frequency"] = str(round(float(A + C) / N, 3))
    snp_2["allele_1"] = snp_2_allele_1

    snp_2_allele_2 = {}
    snp_2_allele_2["allele"] = sorted(hap)[1].split("_")[1]
    snp_2_allele_2["count"] = str(B + D)
    snp_2_allele_2["frequency"] = str(round(float(B + D) / N, 3))
    snp_2["allele_2"] = snp_2_allele_2
    output["snp2"] = snp_2

    two_by_two = {}
    cells = {}
    cells["c11"] = str(A)
    cells["c12"] = str(B)
    cells["c21"] = str(C)
    cells["c22"] = str(D)
    two_by_two["cells"] = cells
    two_by_two["total"] = str(N)
    output["two_by_two"] = two_by_two

    haplotypes = {}
    hap_1 = {}
    hap_1["alleles"] = hap1
    hap_1["count"] = str(hap[hap1])
    hap_1["frequency"] = str(round(float(hap[hap1]) / N, 3))
    haplotypes["hap1"] = hap_1

    hap_2 = {}
    hap_2["alleles"] = hap2
    hap_2["count"] = str(hap[hap2])
    hap_2["frequency"] = str(round(float(hap[hap2]) / N, 3))
    haplotypes["hap2"] = hap_2

    hap_3 = {}
    hap_3["alleles"] = hap3
    hap_3["count"] = str(hap[hap3])
    hap_3["frequency"] = str(round(float(hap[hap3]) / N, 3))
    haplotypes["hap3"] = hap_3

    hap_4 = {}
    hap_4["alleles"] = hap4
    hap_4["count"] = str(hap[hap4])
    hap_4["frequency"] = str(round(float(hap[hap4]) / N, 3))
    haplotypes["hap4"] = hap_4
    output["haplotypes"] = haplotypes

    statistics = {}
    if Ms != 0:
        statistics["d_prime"] = str(round(D_prime, 4))
        statistics["r2"] = str(round(r2, 4))
        statistics["chisq"] = str(round(chisq, 4))
        if p <= 0.00001:
            statistics["p"] = str(round(p, 4))
        else:
            statistics["p"] = "<0.00001"
    else:
        statistics["d_prime"] = D_prime
        statistics["r2"] = r2
        statistics["chisq"] = chisq
        statistics["p"] = p

    output["statistics"] = statistics

    output["corr_alleles"] = corr_alleles
    output["corr_dictionary"] = corr_lookup

    output['populations'] = pops

    # Generate output file
    if write_summary:
        if write_summary == 'STDOUT':
            ldpair_out = sys.stdout
        elif write_summary == 'STDERR':
            ldpair_out = sys.stderr
        else:
            ldpair_out = open(write_summary, 'w')

        print_summary(output, ldpair_out)

        if write_summary != 'STDOUT' and write_summary != 'STDERR':
            ldpair_out.close()

    # Return output
    if return_json:
        return output_json(output)
    else:
        return output

def print_summary(json_dict, outfile=sys.stdout):
    """Write a formatted summary to file from a dictionary.

    Parameters
    ----------
    json_dict : dict
        A dictionary of the kind generated by calculate_pair()
    filehandle : IOBase, optional
        An open writable file handle
    """
    outfile.write("Query SNPs:\n")
    outfile.write(json_dict["snp1"]["rsnum"] + " (" + json_dict["snp1"]["coord"] + ")\n")
    outfile.write(json_dict["snp2"]["rsnum"] + " (" + json_dict["snp2"]["coord"] + ")\n")
    outfile.write("\n")
    outfile.write(','.join(json_dict['populations']) + " Haplotypes:\n")
    outfile.write(" " * 15 + json_dict["snp2"]["rsnum"] + '\n')
    outfile.write(" " * 15 + json_dict["snp2"]["allele_1"]["allele"] + " " * 7 + json_dict["snp2"]["allele_2"]["allele"] + '\n')
    outfile.write(" " * 13 + "-" * 17 + '\n')
    outfile.write(" " * 11 + json_dict["snp1"]["allele_1"]["allele"] + " | " + json_dict["two_by_two"]["cells"]["c11"] + " " * (5 - len(json_dict["two_by_two"]["cells"]["c11"])) + " | " + json_dict["two_by_two"]["cells"]["c12"] + " " * (5 - len(json_dict["two_by_two"]["cells"]["c12"])) + " | " + json_dict["snp1"]["allele_1"]["count"] + " " * (5 - len(json_dict["snp1"]["allele_1"]["count"])) + " (" + json_dict["snp1"]["allele_1"]["frequency"] + ")\n")
    outfile.write(json_dict["snp1"]["rsnum"] + " " * (10 - len(json_dict["snp1"]["rsnum"])) + " " * 3 + "-" * 17 + '\n')
    outfile.write(" " * 11 + json_dict["snp1"]["allele_2"]["allele"] + " | " + json_dict["two_by_two"]["cells"]["c21"] + " " * (5 - len(json_dict["two_by_two"]["cells"]["c21"])) + " | " + json_dict["two_by_two"]["cells"]["c22"] + " " * (5 - len(json_dict["two_by_two"]["cells"]["c22"])) + " | " + json_dict["snp1"]["allele_2"]["count"] + " " * (5 - len(json_dict["snp1"]["allele_2"]["count"])) + " (" + json_dict["snp1"]["allele_2"]["frequency"] + ")\n")
    outfile.write(" " * 13 + "-" * 17 + '\n')
    outfile.write(" " * 15 + json_dict["snp2"]["allele_1"]["count"] + " " * (5 - len(json_dict["snp2"]["allele_1"]["count"])) + " " * 3 + json_dict["snp2"]["allele_2"]["count"] + " " * (5 - len(json_dict["snp2"]["allele_2"]["count"])) + " " * 3 + json_dict["two_by_two"]["total"] + '\n')
    outfile.write(" " * 14 + "(" + json_dict["snp2"]["allele_1"]["frequency"] + ")" + " " * (5 - len(json_dict["snp2"]["allele_1"]["frequency"])) + " (" + json_dict["snp2"]["allele_2"]["frequency"] + ")" + " " * (5 - len(json_dict["snp2"]["allele_2"]["frequency"])) + '\n')
    outfile.write("\n")
    outfile.write("          " + json_dict["haplotypes"]["hap1"]["alleles"] + ": " + json_dict["haplotypes"]["hap1"]["count"] + " (" + json_dict["haplotypes"]["hap1"]["frequency"] + ")\n")
    outfile.write("          " + json_dict["haplotypes"]["hap2"]["alleles"] + ": " + json_dict["haplotypes"]["hap2"]["count"] + " (" + json_dict["haplotypes"]["hap2"]["frequency"] + ")\n")
    outfile.write("          " + json_dict["haplotypes"]["hap3"]["alleles"] + ": " + json_dict["haplotypes"]["hap3"]["count"] + " (" + json_dict["haplotypes"]["hap3"]["frequency"] + ")\n")
    outfile.write("          " + json_dict["haplotypes"]["hap4"]["alleles"] + ": " + json_dict["haplotypes"]["hap4"]["count"] + " (" + json_dict["haplotypes"]["hap4"]["frequency"] + ")\n")
    outfile.write("\n")
    outfile.write("          D': " + json_dict["statistics"]["d_prime"] + '\n')
    outfile.write("          R2: " + json_dict["statistics"]["r2"] + '\n')
    outfile.write("      Chi-sq: " + json_dict["statistics"]["chisq"] + '\n')
    outfile.write("     p-value: " + json_dict["statistics"]["p"] + '\n')
    outfile.write("\n")
    if len(json_dict["corr_alleles"]) == 2:
        outfile.write(json_dict["corr_alleles"][0] + '\n')
        outfile.write(json_dict["corr_alleles"][1] + '\n')
    else:
        outfile.write(json_dict["corr_alleles"][0] + '\n')

    if 'warning' in json_dict:
        outfile.write("WARNING: " + json_dict["warning"] + "!\n")


def get_arg_parser(argv=None):
    """Return command line parser."""

    parser  = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Positional arguments
    parser.add_argument('snp1', help="First SNP (rsID or chr:pos)")
    parser.add_argument('snp2', help="Second SNP (rsID or chr:pos)")
    parser.add_argument('populations',
                        help="Comma separated list of populations")
    parser.add_argument('request', nargs='?',
                        help="String to use for output file name")

    return parser


def main(argv=None):
    """Run as a script."""
    if not argv:
        argv = sys.argv[1:]

    parser = get_arg_parser()
    args = parser.parse_args(argv)

    pops = args.populations.split(',')
    if args.request:
        outfile = "LDpair_" + args.request + ".txt"
    else:
        outfile = None

    # Run function
    json_dict = calculate_pair(args.snp1, args.snp2, pops, outfile)

    # Print output
    if 'error' in json_dict:
        sys.stderr.write('\n\n' + json_dict["error"] + '\n\n')
    else:
        print_summary(json_dict)

if __name__ == '__main__' and '__file__' in globals():
    sys.exit(main())
