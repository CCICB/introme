#!/usr/bin/env Rscript

### Script to generate consensus scores for Introme

if(!require("ROCR")) {
  install.packages("ROCR", repos = "https://cloud.r-project.org/", quiet = TRUE)
  suppressPackageStartupMessages(library(ROCR))
}

if(!require("C50")) {
  install.packages(install.packages("C50", repos = "https://cloud.r-project.org/", dependencies = TRUE))
  suppressPackageStartupMessages(library(C50))
}

if(!require("caret")) {
  install.packages(install.packages("caret", repos = "https://cloud.r-project.org/", dependencies = TRUE))
  suppressPackageStartupMessages(library(caret))
}

# Read in running parameters
args = commandArgs(trailingOnly=TRUE)
Introme_file <- read.csv(args[2], header=T, sep="\t", na.strings = '.')
Introme_model <- readRDS("models/C50_210826.rds")

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
Introme_file$Branchpointer_max_Prob <- as.numeric(Introme_file$Branchpointer_max_Prob)
Introme_file$Branchpointer_U2_Binding_Energy <- as.numeric(Introme_file$Branchpointer_U2_Binding_Energy)
Introme_file$Branchpointer_options <- as.numeric(Introme_file$Branchpointer_options)
Introme_file$Branchpointer_max_U2_Binding_Energy <- as.numeric(Introme_file$Branchpointer_max_U2_Binding_Energy)
Introme_file$Intron_Type <- as.character(Introme_file$Intron_Type)
Introme_file$Intron_Type[is.na(Introme_file$Intron_Type)] <- "U2"
Introme_file$Intron_Type <- as.factor(Introme_file$Intron_Type)
Introme_file$Gene_Regions <- as.factor(Introme_file$Gene_Regions)

# Calculate Introme scores
Predict_scores <- predict(Introme_model, newdata = Introme_file, type = "prob", na.action = na.pass)
Introme_file$Introme <- Predict_scores$SAV

# To toggle: Threshold filtering
# Introme_file_thresholds <- subset(Introme_file, Introme_file$Introme >= 0.52)
# write.table(Introme_file_thresholds, file=args[3]_threshold, quote=FALSE, sep='\t', row.names=FALSE, na = ".")

# Save output
write.table(Introme_file, file=args[3], quote=FALSE, sep='\t', row.names=FALSE, na = ".")
