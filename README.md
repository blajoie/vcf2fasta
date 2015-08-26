# vcf2fasta
reference fasta to diploid fasta via SNPs in VCF
 
## Usage

```

$ python scripts/vcf2fasta.py --help
usage: vcf2fasta.py [-h] -r REF_FASTA -v VCF_FILE -n NAME [--verbose]
                    [--version]

vcf2fasta (diploid)

optional arguments:
  -h, --help            show this help message and exit
  -r REF_FASTA, --ref REF_FASTA
                        input reference file (fasta) (default: None)
  -v VCF_FILE, --vcf VCF_FILE
                        input vcf file (vcf) (default: None)
  -n NAME, --name NAME  sample name (column header) (default: None)
  --verbose             Increase verbosity (specify multiple times for more)
                        (default: None)
  --version             show program's version number and exit

  
## Usage Examples

```

python scripts/vcf2fasta.py -r test/dm3__chr2L.fa.gz -v test/DGRP-057.vcf -n DGRP-057
python scripts/vcf2fasta.py -r test/dm3__chr2L.fa.gz -v test/DGRP-057.vcf -n DGRP-057 --verbose
 
## Change Log

## Bugs and Feedback

For bugs, questions and discussions please use the [Github Issues](https://github.com/blajoie/single-element-stack/issues).

## LICENSE

Licensed under the Apache License, Version 2.0 (the 'License');
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

<http://www.apache.org/licenses/LICENSE-2.0>

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an 'AS IS' BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

