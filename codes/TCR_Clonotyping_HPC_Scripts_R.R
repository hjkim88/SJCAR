###
#   File name : TCR_Clonotyping_HPC_Scripts_R.R
#   Author    : Hyunjin Kim
#   Date      : Apr 27, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Make R scripts for HPC run for TCR clonotyping of the SJCAR19 data
#
#   Instruction
#               1. Source("TCR_Clonotyping_HPC_Scripts_R.R")
#               2. Run the function "make_scripts" - specify the input file paths and the output directory
#               3. The result script files will be generated under the output directory
#
#   Example
#               > source("The_directory_of_TCR_Clonotyping_HPC_Scripts_R.R/TCR_Clonotyping_HPC_Scripts_R.R")
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
  
  ### result script's first line
  script <- "source(\"/home/hkim8/SJCAR19/codes/TCR_Clonotyping_HPC.R\")\n"
  
  ### create R scripts for the HPC run
  for(lib in libs) {
    for(opt in options) {
      for(gap in c(0, 1, 2)) {
        script2 <- paste0(script, "clonotyping_hpc(metadata_path=\"/home/hkim8/SJCAR19/data/metadata_hpc.rda\",
                          lib=\"", lib, "\",
                          option=\"", opt, "\",
                          gap=", gap, ",
                          outputDir=\"/home/hkim8/SJCAR19/results/\")\n")
        ### save the shell script
        write.table(script2, file=paste0(outputDir, lib, "_", opt, "_", gap, ".R"), sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)
        
        ### dos2unix to the script
        system(paste(dos2unixPath, paste0(outputDir, lib, "_", opt, "_", gap, ".R")))
      }
    }
  }
  
}
