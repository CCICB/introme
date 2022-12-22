# Introme

Introme is an _in silico_ splice predictor which evaluates a variant’s likelihood of altering splicing by combining predictions from multiple splice-scoring tools, combined with additional splicing rules, and gene architecture features. Introme can accurately predict the impact of coding and noncoding variants on splicing through investigating for the potential damage, creation or strengthening of splice elements and outperforms all leading tools that we tested.

## Requirements

### Software requirements
- [Docker](https://docs.docker.com/get-docker/)
- [vcfanno](https://github.com/brentp/vcfanno) 
- [spliceai](https://github.com/Illumina/SpliceAI)
- [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
- [bcftools](http://www.htslib.org/download/)
- [samtools](http://www.htslib.org/download/)
- [htslib](http://www.htslib.org/download/)
- [R](https://www.r-project.org/)
- R packages: ROCR, caret
- python3
- python packages: pysam, csv, Bio.Seq, argparse

### Variant Annotation file requirements
Introme requires the following files to be downloaded and placed in the annotations folder.

- CADD v1.3 VCF created using the instructions at: https://github.com/brentp/vcfanno/blob/master/docs/examples/cadd.md
- gnomad.genomes.sites.merged.AF_AC_AN_only.vcf.gz
- [MGRB](https://www.garvan.org.au/research/kinghorn-centre-for-clinical-genomics/research-programs/sydney-genomics-collaborative/mgrb), [Pinese et al, 2020](https://www.nature.com/articles/s41467-019-14079-0)
- [SPIDEX v1.0](https://www.openbioinformatics.org/annovar/spidex_download_form.php)
- [dbscSNV v1.1](http://www.liulab.science/dbscsnv.html)

### Additional file requirements
- A gtf file, ideally containing only protein coding regions to speed up annotations
- A bed file with the regions of interest, you can use either a restricted gene list or a protein coding bed file

## Installation

### WDL Install
The wdl script is labelled Introme.wdl in the [wdl_scripts](https://github.com/CCICB/introme/tree/master/wdl_scripts) folder. These scripts were set up for implementation using Terra. All of the annotation files are required to be in the same folder, and specified as inputs to ensure proper annotation using vcfanno, these requirements will be further documented in the folder.

### Local Install
1. Install the above software requirements and pull Introme. 
2. Download the required annotation files and file requirements and place them in the annotations folder.
3. Update the _.conf_ files with the correct paths (shouldn't be necessary if the same annotation files are used).
4. Pull [MMSplice](https://github.com/gagneurlab/MMSplice_MTSplice) and [Spliceogen](https://github.com/VCCRI/Spliceogen) into the Introme folder from the links provided (ensure the original docker files are not deleted).
5. Build the docker containers for MMSplice and Spliceogen using the code below. If you tag the containers differently, ensure you update the `docker run` section in the _run_introme.sh_ script.

```
cd MMSplice
docker build -t mmsplice .
```
```
cd Spliceogen
docker build -t spliceogen .
```
Note: The MMSplice Docker Container requires more memory than the standard settings for Docker. Upgrade the memory to 10GB to ensure it runs.

6. Run introme using the run command `./run_introme.sh -r $genome.fa -g $gtf -b $bed -v $input -p $prefix`

### Docker Local Install
A more streamlined install of introme for running locally is being developed using Docker. 

    
## Running Introme
Introme can be run using either a local installation, or Docker. 

Furthermore, we have wrapped Introme in Workflow Description Language and implemented using [Terra](https://terra.bio/). We are currently in the process of implementing Introme using [CAVATICA](https://www.cavatica.org), which uses the [SevenBridges Genomics](https://www.sevenbridges.com/) platform.

### Required parameters
- `g` Input GTF file
- `p` Output file prefix
- `r` Reference genome
- `v` Input VCF file

### Optional parameters
- `a` Genome assembly (can be inferred from genome build if in the file name)
- `b` Input BED file (i.e. regions of interest)
- `f` Score all variants ≤ a specified variant allele frequency
- `q` Score all variants regardless of quality score
- `s` Turn off Introme single score check

### Examples

Run Introme with base parameters:
`./run_introme.sh -r $genome.fa -g $gtf -v $input -p $prefix`

Run Introme on a specified gene list (BED format) for variants below 0.1% allele frequency:
`./run_introme.sh -r $genome.fa -g $gtf -v $input -p $prefix -b $bed -f 0.001`

## Interpreting Introme Results

The variant-level scores and supporting information are then fed into the Introme decision tree model to classify the likelihood of a variant altering splicing, which produces an Introme score from 0–1. **We recommend the use of 0.54 as a threshold**, producing a sensitivity of 0.9 and a specificity of 0.95, calculated on the validation dataset. When high specificity is required, a threshold of 0.75 results in a sensitivity of 0.8 and a specificity of 0.97.

We are working on implementing automatic interpretation for the outcome of the splice-altering variant. Until this feature is in place, all of the input scores which make up Introme's final prediction are included in the final .tsv file if further information on the variant prediction is required. 

## Reference Genome versions
Introme currently supports VCF files aligned to the both GRCh37 and GRCh38 reference genomes. Please specify using -a 'hg19/hg38' if your reference genome is not specified in the name of the fasta file.

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

## Support

For additional questions or assistance using Introme, contact psullivan@ccia.org.au (Patricia Sullivan).
