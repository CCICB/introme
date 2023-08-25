
# python3 ESE_ESS_scoring.py AG_check/files/SpliceVarDB.hg38.strand.vcf ESE/working_files/SpliceVarDB.ESE.big.vcf  AG_check/files/hg38_.fa

# python3 ESE/scoring.py AG_check/files/SpliceVarDB.hg38.strand.vcf ESE/working_files/SpliceVarDB.ESE.big.new.vcf  AG_check/files/hg38_.fa

size=large

if test $size = "medium"; then
    time python3 ESE/scoring.py AG_check/files/SpliceVarDB.hg38.strand.vcf ESE/working_files/mid.${size}.strand.vcf  AG_check/files/hg38_.fa
    time python3 ESE/scoring2.py AG_check/files/SpliceVarDB.hg38.strand.vcf ESE/working_files/late.${size}.strand.vcf  AG_check/files/hg38_.fa

    diff ESE/working_files/mid.${size}.strand.vcf ESE/working_files/late.${size}.strand.vcf

elif test $size = "large"; then
    time python3 ESE/scoring.py AG_check/files/SpliceVarDB.hg38.sort.strand.vcf ESE/working_files/mid.${size}.strand.vcf  AG_check/files/hg38_.fa
    time python3 ESE/scoring2.py AG_check/files/SpliceVarDB.hg38.sort.strand.vcf ESE/working_files/late.${size}.strand.vcf  AG_check/files/hg38_.fa     

    diff ESE/working_files/mid.${size}.strand.vcf  ESE/working_files/late.${size}.strand.vcf
fi

