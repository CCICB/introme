# Introme

Introme is an _in silico_ splice predictor which evaluates a variant’s likelihood of altering splicing by combining predictions from multiple splice-scoring tools, combined with additional splicing rules, and gene architecture features. Introme can accurately predict the impact of coding and noncoding variants on splicing through investigating for the potential damage, creation or strengthening of splice elements and outperforms all leading tools that we tested.

## Requirements

We suggest using the dockerised versions of Introme, ideally implemented on.

Introme has also been wrapped in Workflow Description Language and implemented using [Terra](https://terra.bio/). We are currently in the process of implementing Introme using [CAVATICA](https://www.cavatica.org).


Since Introme combines multiple splice-prediction tools, there is a large number of dependencies required to run Introme locally. 

## Running Introme

### Requirements
- `b` Input BED file (i.e. regions of interest)
- `g` Input GTF file
- `p` Output file prefix
- `r` Reference genome
- `v` Input VCF file

### Options
- `f` Score all variants regardless of variant frequency
- `l` List of genes to filter for (.txt file or list accepted) 
- `q` Score all variants regardless of quality score

### Example

`./run_introme.sh -r $genome.fa -g $gtf -b $bed -v $input -p $prefix`


## Interpreting Introme Outputs

The variant-level scores and supporting information are then fed into the Introme decision tree model to classify the likelihood of a variant altering splicing, which produces an Introme score from 0–1. **We recommend the use of 0.54 as a threshold**, producing a sensitivity of 0.9 and a specificity of 0.95, calculated on the validation dataset. When high specificity is required, a threshold of 0.75 results in a sensitivity of 0.8 and a specificity of 0.97.

We are working on implementing automatic interpretation for the outcome of the splice-altering variant. Until this feature is in place, all of the input scores which make up Introme's final prediction are included in the final .tsv file if further information on the variant prediction is required. 


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
