# Introme

Introme is an _in silico_ splice predictor which evaluates a variant’s likelihood of altering splicing by combining predictions from multiple splice-scoring tools, combined with additional splicing rules, and gene architecture features. Introme can accurately predict the impact of coding and noncoding variants on splicing through investigating for the potential damage, creation or strengthening of splice elements and outperforms all leading tools that we tested.

## Requirements

### Software requirements
- vcfanno
- spliceai
- vep
- MMSplice vep plugin and python package
- bedtools
- bcftools
- samtools
- tabix
- R
- R packages: ROCR, caret
- python3
- python packages: pysam, csv, Bio.Seq, argparse, mmsplice

### Variant Annotation file requirements
- CADD v1.3 VCF created using the instructions at: https://github.com/brentp/vcfanno/blob/master/docs/examples/cadd.md
- gnomad.genomes.sites.merged.AF_AC_AN_only.vcf.gz
- [MGRB](https://www.garvan.org.au/research/kinghorn-centre-for-clinical-genomics/research-programs/sydney-genomics-collaborative/mgrb), [Pinese et al, 2020](https://www.nature.com/articles/s41467-019-14079-0)
- CADD v1.3 
- [SPIDEX v1.0](https://www.openbioinformatics.org/annovar/spidex_download_form.php)
- dbscSNV v1.1

## Installation
We suggest using the dockerised versions of Introme, _instructions to be finalised_.
    
## Running Introme
Introme can be run using either a local installation, or Docker. 

Furthermore, we have wrapped Introme in Workflow Description Language and implemented using [Terra](https://terra.bio/). We are currently in the process of implementing Introme using [CAVATICA](https://www.cavatica.org), which uses the [SevenBridges Genomics](https://www.sevenbridges.com/) platform.

### Required parameters
- `b` Input BED file (i.e. regions of interest)
- `g` Input GTF file
- `p` Output file prefix
- `r` Reference genome
- `v` Input VCF file

### Optional parameters
- `f` Score all variants regardless of variant frequency
- `l` List of genes to filter for (.txt file or list accepted) 
- `q` Score all variants regardless of quality score

### Example

`./run_introme.sh -r $genome.fa -g $gtf -b $bed -v $input -p $prefix`


## Interpreting Introme Results

The variant-level scores and supporting information are then fed into the Introme decision tree model to classify the likelihood of a variant altering splicing, which produces an Introme score from 0–1. **We recommend the use of 0.54 as a threshold**, producing a sensitivity of 0.9 and a specificity of 0.95, calculated on the validation dataset. When high specificity is required, a threshold of 0.75 results in a sensitivity of 0.8 and a specificity of 0.97.

We are working on implementing automatic interpretation for the outcome of the splice-altering variant. Until this feature is in place, all of the input scores which make up Introme's final prediction are included in the final .tsv file if further information on the variant prediction is required. 

## Reference Genome versions
Introme currently supports VCF files aligned to the hs37d5 reference genome. GRCh38 is not yet supported.

## Funding

The development of Introme has been supported grants, fellowships and scholarships provided by:
- Luminesce Alliance
- Cancer Australia and My Room
- NHMRC
- NSW Health
- Australian Government Research Training Program
- The Kids Cancer Alliance
- Petre Foundation
- Fulbright Future Scholarship

## Development

Introme was initially developed by Dr. Mark Cowley, Dr. Velimir Gayevskiy and Dr. Sarah Beecroft at the Garvan Institute's Kinghorn Centre for Clinical Genomics, and the initial implementation can be found at [KCCG's Introme Repository](https://github.com/KCCG/introme). 

Introme has since been adapted and reimplemented by Patricia Sullivan, Dr. Mark Cowley and Dr. Mark Pinese at the Children's Cancer Institute. This version extends on KCCG's Introme in terms of accuracy, the addition of mulitple splice-scoring tools, and the use of machine learning.
