# Import
import argparse
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_save, predict_all_table
from mmsplice.utils import max_varEff
from mmsplice.utils import writeVCF

## Bring in arguments
parser = argparse.ArgumentParser()
parser.add_argument("--vcf", action="store", dest="vcf", help="input gzipped vcf")
parser.add_argument("--fasta", action="store", dest="fasta", help="reference genome fasta file")
parser.add_argument("--gtf", action="store", dest="gtf", help="gtf file")
parser.add_argument("--output", action="store", dest="output", help="output file for MMSplice results")
args = parser.parse_args()

# Specify model
model = MMSplice()

#dl = SplicingVCFDataloader(gtf, fasta, vcf, encode=False, tissue_specific=False)
dl = SplicingVCFDataloader(args.gtf, args.fasta, args.vcf)

# Or predict and return as df
predictions = predict_all_table(model, dl, pathogenicity=True, splicing_efficiency=True)

# Summerize with maximum effect size
predictionsMax = max_varEff(predictions)

writeVCF(args.vcf, args.output, predictionsMax)
