# vcf2fasta
reference fasta to diploid fasta via SNPs in VCF

options
-r (input reference fasta file)
all chromosomes can be concatenated together into single file
gzipped fasta is OK

-v 
VCF file

-n 
column header for the genotype information

--verbose
turn on verbose

usage:
python scripts/vcf2fasta.py -r test/dm3__chr2L.fa.gz -v test/DGRP-057.vcf -n DGRP-057
python scripts/vcf2fasta.py -r test/dm3__chr2L.fa.gz -v test/DGRP-057.vcf -n DGRP-057 --verbose
 
  