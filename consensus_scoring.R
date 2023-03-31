#!/usr/bin/env Rscript

### Script to generate consensus scores for Introme

if(!require("ROCR")) {
  install.packages("ROCR", repos = "http://cran.us.r-project.org", quiet = TRUE)
  suppressPackageStartupMessages(library(ROCR))
}

if(!require("caret")) {
  install.packages("caret", quiet = TRUE)
  suppressPackageStartupMessages(library(caret))
}

# Read in running parameters
args = commandArgs(trailingOnly=TRUE)
Introme_file <- read.csv(args[2], header=T, sep="\t", na.strings = '.')
Introme_model <- readRDS("models/C50_210826.rds")

Introme_file$AG_Created <- NA
Introme_file$GT_Created <- NA
Introme_file$AG_Lost <- NA
Introme_file$GT_Lost <- NA

# Ensure all columns are in the right format (for when there are missing values)
Introme_file$SPIDEX_dPSI_Zscore <- as.numeric(Introme_file$SPIDEX_dPSI_Zscore)
Introme_file$SPIDEX_dPSI_Max_Tissue <- as.numeric(Introme_file$SPIDEX_dPSI_Max_Tissue)
Introme_file$dbscSNV_AdaBoost_Score <- as.numeric(Introme_file$dbscSNV_AdaBoost_Score)
Introme_file$dbscSNV_RandomForest_Score <- as.numeric(Introme_file$dbscSNV_RandomForest_Score)
Introme_file$AG_Created <- as.factor(Introme_file$AG_Created)
Introme_file$GT_Created <- as.factor(Introme_file$GT_Created)
Introme_file$AG_Lost <- as.factor(Introme_file$AG_Lost)
Introme_file$GT_Lost <- as.factor(Introme_file$GT_Lost)
Introme_file$Branchpointer_Prob <- as.numeric(Introme_file$Branchpointer_Prob)
Introme_file$Branchpointer_U2_Binding_Energy <- as.numeric(Introme_file$Branchpointer_U2_Binding_Energy)
Introme_file$Branchpointer_max_Prob <- as.numeric(Introme_file$Branchpointer_max_Prob)
Introme_file$Branchpointer_options <- as.numeric(Introme_file$Branchpointer_options)
Introme_file$Branchpointer_max_U2_Binding_Energy <- as.numeric(Introme_file$Branchpointer_max_U2_Binding_Energy)

Introme_file$Intron_Type[is.na(Introme_file$Intron_Type)] <- "U2"
Introme_file$Intron_Type <- as.factor(Introme_file$Intron_Type)
Introme_file$Gene_Regions <- as.factor(Introme_file$Gene_Regions)
Introme_file$MMSplice_alt_acceptor <- as.numeric(Introme_file$MMSplice_alt_acceptor)
Introme_file$MMSplice_alt_acceptor_intron <- as.numeric(Introme_file$MMSplice_alt_acceptor_intron)
Introme_file$MMSplice_alt_donor <- as.numeric(Introme_file$MMSplice_alt_donor)
Introme_file$MMSplice_alt_donor_intron <- as.numeric(Introme_file$MMSplice_alt_donor_intron)
Introme_file$MMSplice_alt_exon <- as.numeric(Introme_file$MMSplice_alt_exon)
Introme_file$MMSplice_delta_logit_PSI <- as.numeric(Introme_file$MMSplice_delta_logit_PSI)
Introme_file$MMSplice_pathogenicity <- as.numeric(Introme_file$MMSplice_pathogenicity)
Introme_file$MMSplice_ref_acceptor <- as.numeric(Introme_file$MMSplice_ref_acceptor)
Introme_file$MMSplice_ref_acceptor_intron <- as.numeric(Introme_file$MMSplice_ref_acceptor_intron)
Introme_file$MMSplice_ref_donor <- as.numeric(Introme_file$MMSplice_ref_donor)
Introme_file$MMSplice_ref_donor_intron <- as.numeric(Introme_file$MMSplice_ref_donor_intron)
Introme_file$MMSplice_ref_exon <- as.numeric(Introme_file$MMSplice_ref_exon)

# Calculate Introme scores
Predict_scores <- predict(Introme_model, newdata = Introme_file, type = "prob", na.action = na.pass)
Introme_file$Introme <- Predict_scores$SAV

# Remove false positives caused by no scores (caused by lack of non-scores in the negative training data set)
Introme_file[apply(Introme_file[17:56], 1, function(x) all(is.na(x)) == "TRUE"),"Introme"] <- 0

# To toggle: Threshold filtering
# Introme_file_thresholds <- subset(Introme_file, Introme_file$Introme >= 0.61)
# write.table(Introme_file_thresholds, file=args[3]_threshold, quote=FALSE, sep='\t', row.names=FALSE, na = ".")

# Save output
write.table(Introme_file, file=args[3], quote=FALSE, sep='\t', row.names=FALSE, na = ".")
