#!/usr/bin/env Rscript

### Script to generate consensus scores for Introme

require("ROCR")
require("caret")
#install.packages("ROCR", repos = "http://cran.us.r-project.org", quiet = TRUE)
#suppressPackageStartupMessages(library(ROCR))

# Read in running parameters
args = commandArgs(trailingOnly=TRUE)

Introme_file <- read.csv(args[2], header=T, sep="\t", na.strings = '.')

if (args[1] == "full") {
  Introme_model <- readRDS("models/C50_200812.rds")
} else if (args[1] == "fast") {
  Introme_model <- readRDS("models/C50_Fast.rds")
  
  # Assign all entries which were not scored in the SpliceAI Precomputed file with a 0
  Introme_file$SpliceAI_Acceptor_Gain[is.na(Introme_file$SpliceAI_Acceptor_Gain)] <- 0
  Introme_file$SpliceAI_Acceptor_Loss[is.na(Introme_file$SpliceAI_Acceptor_Loss)] <- 0
  Introme_file$SpliceAI_Donor_Gain[is.na(Introme_file$SpliceAI_Donor_Gain)] <- 0
  Introme_file$SpliceAI_Donor_Loss[is.na(Introme_file$SpliceAI_Donor_Loss)] <- 0
}

Introme_file$SPIDEX_dPSI_Zscore <- as.numeric(Introme_file$SPIDEX_dPSI_Zscore)
Introme_file$SPIDEX_dPSI_Max_Tissue <- as.numeric(Introme_file$SPIDEX_dPSI_Max_Tissue)
Introme_file$dbscSNV_AdaBoost_Score <- as.numeric(Introme_file$dbscSNV_AdaBoost_Score)
Introme_file$dbscSNV_RandomForest_Score <- as.numeric(Introme_file$dbscSNV_RandomForest_Score)
Introme_file$AG_Created <- as.numeric(Introme_file$AG_Created)
Introme_file$GT_Created <- as.numeric(Introme_file$GT_Created)
Introme_file$Branchpointer_Branchpoint_Prob <- as.numeric(Introme_file$Branchpointer_Branchpoint_Prob)
Introme_file$Branchpointer_U2_Binding_Energy <- as.numeric(Introme_file$Branchpointer_U2_Binding_Energy)


Predict_scores <- predict(Introme_model, newdata = Introme_file, type = "prob", na.action = na.pass)
Introme_file$Introme <- Predict_scores$SAV

# To toggle: Threshold filtering
# Introme_file_thresholds <- subset(Introme_file, Introme_file$Introme >= 0.52)
# write.table(Introme_file_thresholds, file=args[3]_threshold, quote=FALSE, sep='\t', row.names=FALSE, na = ".")

# Save output
write.table(Introme_file, file=args[3], quote=FALSE, sep='\t', row.names=FALSE, na = ".")
