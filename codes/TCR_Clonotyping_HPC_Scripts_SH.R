###
#   File name : TCR_Clonotyping_HPC_Scripts_SH.R
#   Author    : Hyunjin Kim
#   Date      : Apr 27, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Make SHELL scripts for HPC run for TCR clonotyping of the SJCAR19 data
#
#   Instruction
#               1. Source("TCR_Clonotyping_HPC_Scripts_SH.R")
#               2. Run the function "make_scripts" - specify the input file paths and the output directory
#               3. The result script files will be generated under the output directory
#
#   Example
#               > source("The_directory_of_TCR_Clonotyping_HPC_Scripts_SH.R/TCR_Clonotyping_HPC_Scripts_SH.R")
#               > make_scripts(metadata_path="./data/metadata_hpc.rda",
#                              dos2unixPath="C:/Program Files/dos2unix-7.4.1-win64/bin/dos2unix.exe",
#                              outputDir="./codes/hpc_scripts/")
###

make_scripts <- function(metadata_path="./data/metadata_hpc.rda",
                         dos2unixPath="dos2unix.exe",
                         outputDir="./codes/hpc_scripts/") {
  
  ### load the data
  load(metadata_path)
  
  ### unique libraries
  libs <- unique(metadata$Library)
  
  ### options
  options <- c("ab_strict", "a_strict", "b_strict", "ab_lenient")
  
  ### master script's first two lines
  masterScript <- "#!/bin/bash\n"
  
  ### create R scripts for the HPC run
  for(lib in libs) {
    for(opt in options) {
      for(gap in c(0, 1, 2)) {
        ### result script's first line
        script <-"#!/bin/bash\n"
        
        ### complete the shell script
        script <- paste0(script, "export R_LIBS=/home/hkim8/R/x86_64-pc-linux-gnu-library/3.6\n")
        script <- paste0(script, "module load R/3.6.3\n")
        script <- paste0(script, "R_FILE=/home/hkim8/SJCAR19/codes/", lib, "_", opt, "_", gap, ".R\n")
        script <- paste0(script, "R CMD BATCH ${R_FILE}\n")
        
        ### append the master script
        masterScript <- paste0(masterScript, "bsub -R \"rusage[mem=1000]\" -P SJCAR19 -q standard -J ",
                               lib, "_", opt, "_", gap, " \"/home/hkim8/SJCAR19/codes/",
                               lib, "_", opt, "_", gap, ".sh\"\n")
        
        ### save the shell script
        write.table(script, file=paste0(outputDir, lib, "_", opt, "_", gap, ".sh"), sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)
        
        ### dos2unix to the script
        system(paste(dos2unixPath, paste0(outputDir, lib, "_", opt, "_", gap, ".sh")))
      }
    }
  }
  
  ### save the master script
  write.table(masterScript, file=paste0(outputDir, "master.sh"), sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  ### dos2unix to the master script
  system(paste(dos2unixPath, paste0(outputDir, "master.sh")))
  
}
