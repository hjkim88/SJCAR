#!/bin/bash

R_FILE=/research/sharedresources/immunoinformatics/hkim8/SJCAR19_classifier/102821_SJCAR19_Classifier_HPC_R_Script.R

export R_LIBS=/home/hkim8/R/x86_64-pc-linux-gnu-library/4.1

module load R/4.1.0

R CMD BATCH ${R_FILE}
