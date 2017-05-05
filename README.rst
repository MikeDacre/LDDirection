############
LD_Direction
############

Get information about direction of agreement for any two lists of SNPs, based on `LDlink <https://github.com/CBIIT/nci-webtools-dceg-linkage>`_ and `plink <https://www.cog-genomics.org/plink>`_

Compares every SNP in one list to each SNP in another:

.. image:: LD_Direction_graphic.png

Pairs are the filtered by distance on the genome and linkage disequilibrium. The filtered list includes phasing information, allowing instant lookup of the allele in one SNP given an allele in the other. Lookup for a single pair is also possible.

The software works best used as a python library, but two scripts are provided also:

- ``ldlists`` :  allows comparison of two SNP lists as described above
- ``ldpair``  :  allows comparison of just one pair with a more detailed report than given by ``ldlists``

To run this code, you will need a copy of a dbSNP sqlite database and a copy of all the 1000genomes VCF files.

**Note**: Currently under development and unstable/unusable.

The code can be used as 

...........
Requirments
...........

Code
====

To run ``ldlists`` you will need a funtional version of PLINK 1.9: https://www.cog-genomics.org/plink, to run ``ldpair`` you will need a functional version of tabix `htslib <https://github.com/samtools/htslib>`_

I recommend installing both of these tools.

Python Requirements
-------------------

All python requirements will be automatically downloaded when installing using the `setup.py` file.

Essential requirements:

- `Cython <http://cython.org/>`_
- `dbSNP <https://github.com/MikeDacre/dbSNP>`_

Required for table handling in `ldlists`:

- `pandas <http://pandas.pydata.org/>`_

Recommended for progress bars:

- `tqdm <https://pypi.python.org/pypi/tqdm>`_

Recommended if you have access to a torque or slurm cluster:

- `fyrd <https://fyrd.science>`_

All of the above except fyrd will be installed automatically, if you want to parallelize on a cluster, install and configure fyrd on your system.

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
