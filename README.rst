############
LD_Direction
############

Get information about direction of agreement for any two lists of SNPs, based on https://github.com/CBIIT/nci-webtools-dceg-linkage

This basically just wraps the LDpair.py tool from that repository, all rights are maintained by the original developers.

To run this code, you will need a copy of a dbSNP sqlite database and a copy of all the 1000genomes VCF files.

**Note**: Currently under development and unstable/unusable.

...........
Requirments
...........

Code
====

You will need a funtional version of PLINK 1.9: https://www.cog-genomics.org/plink

Python Requirements
-------------------

All python requirements will be automatically downloaded when installing using the `setup.py` file

Data
====

To run this code you will need the 1000genomes data in plink format and a copy
of dbSNP. It could hypothetically be generalized to run on any plink file and
without dbSNP.

1000genomes
-----------

To get the 1000genomes data, download all the `ALL.chr*vcf.gz` files on
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502.

To convert to plink, you will need to run the following commands on every file:

.. code:: shell

  i=ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz;
  j=$(echo $i | sed 's/.vcf.gz//');
  vftools --gzvcf $i --plink-tped --out $j &&
  plink --tped $j.tped --tfam $j.tfam --out $j

To run these on every file:

.. code:: shell

  for i in ALL.chr*.vcf.gz; do
    j=$(echo $i | sed 's/.vcf.gz//');
    vftools --gzvcf $i --plink-tped --out $j &&
    plink --tped $j.tped --tfam $j.tfam --out $j
  done
 
dbSNP
-----

We use a dbSNP database interface that I wrote for position lookup:
https://github.com/MikeDacre/dbSNP You can download any version of dbSNP you
want and build a database using that code. However, this is a very slow
process. You can find pre-built databases at sftp
dbsnp@esme.stanford.edu/dbSNP, password is 'givemethedb'. Note that we support
sftp only, not regular FTP, and this is just my personal machine at Stanford,
if there are too many downloads I will have to stop providing data.
