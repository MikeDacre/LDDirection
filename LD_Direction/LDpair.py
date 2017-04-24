#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Modified LDpair from:

https://github.com/CBIIT/nci-webtools-dceg-linkage/blob/master/LDlink/LDpair.py
"""
import os
import json
import math
import time
import subprocess
from . import dbsnp_db

DB_SNP_VER = '149'

POPULATIONS = ["ALL", "AFR", "AMR", "EAS", "EUR", "SAS", "ACB", "ASW", "BEB",
               "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH",
               "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL",
               "PJL", "PUR", "STU", "TSI", "YRI"]

# Set data directories
DATA_DIR = "/local/content/ldlink/data/"
SNP_DIR = DATA_DIR + "snp142/snp142_annot_2.db"
POP_DIR = DATA_DIR + "1000G/Phase3/samples/"
VCF_DIR = DATA_DIR + "1000G/Phase3/genotypes/ALL.chr"
TMP_DIR = "./tmp/"

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
    """Use dbSNP and 1000genomes tabix to get SNP information."""
    # Find RS numbers in dbSNP
    db = dbsnp_db.DB(version=DB_SNP_VER, location=DATA_DIR)
    snp_coord = db.lookup_rsids(snp)

    # Get SNP info from 1000genomes
    vcf_file = (
        VCF_DIR + snp_coord.chrom +
        ".phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz"
    )
    tabix_cmd = "tabix {0} {1}:{2}-{2} | grep -v -e END".format(
        vcf_file, snp_coord.chrom, snp_coord.start
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

    return snp_coord, geno, allele, head, snp_a1, snp_a2


def calculate_pair(snp1, snp2, pops, outfile):
    """Find LD information for any two SNPs.

    Parameters
    ----------
    snp1/snp2 : str
        rsIDs to compare
    pops : list
        list of populations (e.g. YRI, ESN).
    outfile : str
        a file to write output to

    Returns
    -------
    SNP_Pair
        A SNP_Pair class with information about SNP linkage.
    """
    # Ensure tmp directory exists
    if not os.path.exists(TMP_DIR):
        os.makedirs(TMP_DIR)

    # Create JSON output
    output = {}

    # Get SNP info
    try:
        (snp1_coord, geno1,
         allele1, head1, snp1_a1, snp1_a2) = get_snp_info(snp1)
        (snp2_coord, geno2,
         allele2, head2, snp2_a1, snp2_a2) = get_snp_info(snp2)
    except SNP_Lookup_Failure as e:
        output["error"] = str(e)
        return output_json(output)

    for snp_coord in [snp1_coord, snp2_coord]:
        if snp_coord is None:
            output["error"] = '{} is not in dbSNP build {}.'.format(
                snp_coord.name, DB_SNP_VER
            )
            return output_json(output)

    # Check if SNPs are on the same chromosome
    if snp1_coord.chrom != snp2_coord.chrom:
        output["warning"] = ("{} and {} are on different chromosomes"
                             .format(snp1, snp2))

    for geno, snp in [(geno1, snp1), (geno2, snp2)]:
        if geno[2] != snp:
            if "warning" not in output:
                output["warning"] = ''
            else:
                output["warning"] += '. '
            output["warning"] += (
                "Genomic position for query variant1 ({}) does not match RS "
                "number at 1000G position ({})"
            ).format(snp, geno[2])

    if geno1[1] != snp1_coord[2]:
        output["error"] = "VCF File does not match variant coordinates for SNP1."
        return output_json(output)
    if geno2[1] != snp2_coord[2]:
        output["error"] = "VCF File does not match variant coordinates for SNP2."
        return output_json(output)

    ### Combine phased genotypes ###
    geno = {}
    for i in range(9, len(head1)):
        geno[head1[i]] = [allele1[geno1[i]], ".."]

    for i in range(9, len(head2)):
        if head2[i] in geno:
            geno[head2[i]][1] = allele2[geno2[i]]

    # Select desired ancestral populations
    pop_dirs = []
    bad_pops = []
    for pop_i in pops:
        if pop_i in POPULATIONS:
            pop_dirs.append(POP_DIR + pop_i + ".txt")
        else:
            bad_pops.append(pop_i)
    if bad_pops:
        output["error"] = (
            "{} are not ancestral populations. Choose one of the following "
            "ancestral populations: {}"
        ).format(bad_pops, POPULATIONS)
        return output_json(output)

    get_pops_cmd = "cat {}".format(" ".join(pop_dirs))
    pop_ids = list(set(run_cmnd(get_pops_cmd)))

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
        if dA == dB == dC == dD:
            if tmax == A or tmax == D:
                corr1 = outstr.format(
                    snp1, sorted(hap)[0].split("_")[0],
                    snp2, sorted(hap)[0].split("_")[1]
                )
                corr2 = outstr.format(
                    snp1, sorted(hap)[2].split("_")[0],
                    snp2, sorted(hap)[1].split("_")[1]
                )
                corr_alleles = [corr1, corr2]
            else:
                corr1 = outstr.format(
                    snp1, sorted(hap)[0].split("_")[0],
                    snp2, sorted(hap)[1].split("_")[1]
                )
                corr2 = outstr.format(
                    snp1, sorted(hap)[2].split("_")[0],
                    snp2, sorted(hap)[0].split("_")[1]
                )
                corr_alleles = [corr1, corr2]
        elif dmax == dA or dmax == dD:
            corr1 = outstr.format(
                snp1, sorted(hap)[0].split("_")[0],
                snp2, sorted(hap)[0].split("_")[1]
            )
            corr2 = outstr.format(
                snp1, sorted(hap)[2].split("_")[0],
                snp2, sorted(hap)[1].split("_")[1]
            )
            corr_alleles = [corr1, corr2]
        else:
            corr1 = outstr.format(
                snp1, sorted(hap)[0].split("_")[0],
                snp2, sorted(hap)[1].split("_")[1]
            )
            corr2 = outstr.format(
                snp1, sorted(hap)[2].split("_")[0],
                snp2, sorted(hap)[0].split("_")[1]
            )
            corr_alleles = [corr1, corr2]
    else:
        corr_alleles = [snp1 + " and " + snp2 + " are in linkage equilibrium"]

    # Create JSON output
    snp_1 = {}
    snp_1["rsnum"] = snp1
    snp_1["coord"] = "chr" + snp1_coord[1] + ":" + snp1_coord[2]

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
    snp_2["rsnum"] = snp2
    snp_2["coord"] = "chr" + snp2_coord[1] + ":" + snp2_coord[2]

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
        if p >= 0.0001:
            statistics["p"] = str(round(p, 4))
        else:
            statistics["p"] = "<0.0001"
    else:
        statistics["d_prime"] = D_prime
        statistics["r2"] = r2
        statistics["chisq"] = chisq
        statistics["p"] = p

    output["statistics"] = statistics

    output["corr_alleles"] = corr_alleles

    # Generate output file
    with open(outfile, 'w') as ldpair_out:
        print("Query SNPs:", file=ldpair_out)
        print(output["snp1"]["rsnum"] + " (" + output["snp1"]["coord"] + ")", file=ldpair_out)
        print(output["snp2"]["rsnum"] + " (" + output["snp2"]["coord"] + ")", file=ldpair_out)
        print("", file=ldpair_out)
        print(pops + " Haplotypes:", file=ldpair_out)
        print(" " * 15 + output["snp2"]["rsnum"], file=ldpair_out)
        print(" " * 15 + output["snp2"]["allele_1"]["allele"] + " " * 7 + output["snp2"]["allele_2"]["allele"], file=ldpair_out)
        print(" " * 13 + "-" * 17, file=ldpair_out)
        print(" " * 11 + output["snp1"]["allele_1"]["allele"] + " | " + output["two_by_two"]["cells"]["c11"] + " " * (5 - len(output["two_by_two"]["cells"]["c11"])) + " | " + output["two_by_two"]["cells"]["c12"] + " " * (5 - len(output["two_by_two"]["cells"]["c12"])) + " | " + output["snp1"]["allele_1"]["count"] + " " * (5 - len(output["snp1"]["allele_1"]["count"])) + " (" + output["snp1"]["allele_1"]["frequency"] + ")", file=ldpair_out)
        print(output["snp1"]["rsnum"] + " " * (10 - len(output["snp1"]["rsnum"])) + " " * 3 + "-" * 17, file=ldpair_out)
        print(" " * 11 + output["snp1"]["allele_2"]["allele"] + " | " + output["two_by_two"]["cells"]["c21"] + " " * (5 - len(output["two_by_two"]["cells"]["c21"])) + " | " + output["two_by_two"]["cells"]["c22"] + " " * (5 - len(output["two_by_two"]["cells"]["c22"])) + " | " + output["snp1"]["allele_2"]["count"] + " " * (5 - len(output["snp1"]["allele_2"]["count"])) + " (" + output["snp1"]["allele_2"]["frequency"] + ")", file=ldpair_out)
        print(" " * 13 + "-" * 17, file=ldpair_out)
        print(" " * 15 + output["snp2"]["allele_1"]["count"] + " " * (5 - len(output["snp2"]["allele_1"]["count"])) + " " * 3 + output["snp2"]["allele_2"]["count"] + " " * (5 - len(output["snp2"]["allele_2"]["count"])) + " " * 3 + output["two_by_two"]["total"], file=ldpair_out)
        print(" " * 14 + "(" + output["snp2"]["allele_1"]["frequency"] + ")" + " " * (5 - len(output["snp2"]["allele_1"]["frequency"])) + " (" + output["snp2"]["allele_2"]["frequency"] + ")" + " " * (5 - len(output["snp2"]["allele_2"]["frequency"])), file=ldpair_out)
        print("", file=ldpair_out)
        print("          " + output["haplotypes"]["hap1"]["alleles"] + ": " + output["haplotypes"]["hap1"]["count"] + " (" + output["haplotypes"]["hap1"]["frequency"] + ")", file=ldpair_out)
        print("          " + output["haplotypes"]["hap2"]["alleles"] + ": " + output["haplotypes"]["hap2"]["count"] + " (" + output["haplotypes"]["hap2"]["frequency"] + ")", file=ldpair_out)
        print("          " + output["haplotypes"]["hap3"]["alleles"] + ": " + output["haplotypes"]["hap3"]["count"] + " (" + output["haplotypes"]["hap3"]["frequency"] + ")", file=ldpair_out)
        print("          " + output["haplotypes"]["hap4"]["alleles"] + ": " + output["haplotypes"]["hap4"]["count"] + " (" + output["haplotypes"]["hap4"]["frequency"] + ")", file=ldpair_out)
        print("", file=ldpair_out)
        print("          D': " + output["statistics"]["d_prime"], file=ldpair_out)
        print("          R2: " + output["statistics"]["r2"], file=ldpair_out)
        print("      Chi-sq: " + output["statistics"]["chisq"], file=ldpair_out)
        print("     p-value: " + output["statistics"]["p"], file=ldpair_out)
        print("", file=ldpair_out)
        if len(output["corr_alleles"]) == 2:
            print(output["corr_alleles"][0], file=ldpair_out)
            print(output["corr_alleles"][1], file=ldpair_out)
        else:
            print(output["corr_alleles"][0], file=ldpair_out)

        try:
            output["warning"]
        except KeyError:
            www="do nothing"
        else:
            print("WARNING: " + output["warning"] + "!", file=ldpair_out)

    # Return output
    return(json.dumps(output, sort_keys=True, indent=2))


def main():
    import json
    import sys

    # Import LDpair options
    if len(sys.argv) == 5:
        snp1 = sys.argv[1]
        snp2 = sys.argv[2]
        pop = sys.argv[3]
        request = sys.argv[4]
    elif sys.argv[4] is False:
        snp1 = sys.argv[1]
        snp2 = sys.argv[2]
        pop = sys.argv[3]
        request = str(time.strftime("%I%M%S"))
    else:
        print("Correct useage is: LDpair.py snp1 snp2 populations request")
        sys.exit()

    # Run function
    out_json = calculate_pair(snp1, snp2, pop, request)

    # Print output
    json_dict = json.loads(out_json)
    try:
        json_dict["error"]

    except KeyError:
        print("")
        print("Query SNPs:")
        print(json_dict["snp1"]["rsnum"] + " (" + json_dict["snp1"]["coord"] + ")")
        print(json_dict["snp2"]["rsnum"] + " (" + json_dict["snp2"]["coord"] + ")")
        print("")
        print(pop + " Haplotypes:")
        print(" " * 15 + json_dict["snp2"]["rsnum"])
        print(" " * 15 + json_dict["snp2"]["allele_1"]["allele"] + " " * 7 + json_dict["snp2"]["allele_2"]["allele"])
        print(" " * 13 + "-" * 17)
        print(" " * 11 + json_dict["snp1"]["allele_1"]["allele"] + " | " + json_dict["two_by_two"]["cells"]["c11"] + " " * (5 - len(json_dict["two_by_two"]["cells"]["c11"])) + " | " + json_dict["two_by_two"]["cells"]["c12"] + " " * (5 - len(json_dict["two_by_two"]["cells"]["c12"])) + " | " + json_dict["snp1"]["allele_1"]["count"] + " " * (5 - len(json_dict["snp1"]["allele_1"]["count"])) + " (" + json_dict["snp1"]["allele_1"]["frequency"] + ")")
        print(json_dict["snp1"]["rsnum"] + " " * (10 - len(json_dict["snp1"]["rsnum"])) + " " * 3 + "-" * 17)
        print(" " * 11 + json_dict["snp1"]["allele_2"]["allele"] + " | " + json_dict["two_by_two"]["cells"]["c21"] + " " * (5 - len(json_dict["two_by_two"]["cells"]["c21"])) + " | " + json_dict["two_by_two"]["cells"]["c22"] + " " * (5 - len(json_dict["two_by_two"]["cells"]["c22"])) + " | " + json_dict["snp1"]["allele_2"]["count"] + " " * (5 - len(json_dict["snp1"]["allele_2"]["count"])) + " (" + json_dict["snp1"]["allele_2"]["frequency"] + ")")
        print(" " * 13 + "-" * 17)
        print(" " * 15 + json_dict["snp2"]["allele_1"]["count"] + " " * (5 - len(json_dict["snp2"]["allele_1"]["count"])) + " " * 3 + json_dict["snp2"]["allele_2"]["count"] + " " * (5 - len(json_dict["snp2"]["allele_2"]["count"])) + " " * 3 + json_dict["two_by_two"]["total"])
        print(" " * 14 + "(" + json_dict["snp2"]["allele_1"]["frequency"] + ")" + " " * (5 - len(json_dict["snp2"]["allele_1"]["frequency"])) + " (" + json_dict["snp2"]["allele_2"]["frequency"] + ")" + " " * (5 - len(json_dict["snp2"]["allele_2"]["frequency"])))
        print("")
        print("          " + json_dict["haplotypes"]["hap1"]["alleles"] + ": " + json_dict["haplotypes"]["hap1"]["count"] + " (" + json_dict["haplotypes"]["hap1"]["frequency"] + ")")
        print("          " + json_dict["haplotypes"]["hap2"]["alleles"] + ": " + json_dict["haplotypes"]["hap2"]["count"] + " (" + json_dict["haplotypes"]["hap2"]["frequency"] + ")")
        print("          " + json_dict["haplotypes"]["hap3"]["alleles"] + ": " + json_dict["haplotypes"]["hap3"]["count"] + " (" + json_dict["haplotypes"]["hap3"]["frequency"] + ")")
        print("          " + json_dict["haplotypes"]["hap4"]["alleles"] + ": " + json_dict["haplotypes"]["hap4"]["count"] + " (" + json_dict["haplotypes"]["hap4"]["frequency"] + ")")
        print("")
        print("          D': " + json_dict["statistics"]["d_prime"])
        print("          R2: " + json_dict["statistics"]["r2"])
        print("      Chi-sq: " + json_dict["statistics"]["chisq"])
        print("     p-value: " + json_dict["statistics"]["p"])
        print("")
        if len(json_dict["corr_alleles"]) == 2:
            print(json_dict["corr_alleles"][0])
            print(json_dict["corr_alleles"][1])
        else:
            print(json_dict["corr_alleles"][0])

        try:
            json_dict["warning"]
        except KeyError:
            print("")
        else:
            print("WARNING: " + json_dict["warning"] + "!")
            print("")

    else:
        print("")
        print(json_dict["error"])
        print("")


if __name__ == "__main__":
    main()
