##################################################################################
##                                                                              ##
##                                                                              ##
##                  SCRIPT A: BASECALLING, DEMULTIPLEXING AND                   ##
##                 MERGING FASTQ FILES IN A SINGLE FASTA FILE                   ##
##                                                                              ##
##                                                                              ##
##################################################################################





##### PART 1 - Local computer (bash): Merging all pod5 together no matter if they passed the Q-score filter or not #####

# Note: MinKNOW settings were initially set to a fast basecalling using GUPPY, with a minimum Q-score of 8
# We merged all POD5 (different folders) not to filter per Q, and used the latest sup basecalling model in Dorado

# Note 2: The pod5 were all grouped in the "merged_pod5" so they can be merged, but this folder as well as the "MinKNOW_output" were removed after (too heavy)

# sudo umount /mnt/h
# sudo mount -t drvfs H: /mnt/h
# cd mnt/h/seabass_edna_methylation_data/A_data
# pod5 merge MinKNOW_output/no_sample/20230629_0910_MN24659_FAX01677_cb75b57b/pod5/pod5_all -o merged_pod5/unfiltered_seabass_edna_reads.pod5 --duplicate-ok





##### PART 2 - GPU server (docker + bash): Superior basecalling with Dorado on the full POD5 #####

# Info: Docker, and Dorado (as well as the basecalling model used) were all previously installer on the GPU server

# Note 3: As for the POD5 file, the FASTQ file was not saved because it was too heavy (folder basecalling_sup removed)

# docker run -it -v /home/elruiz/:/home/elruiz/ -v /etc/passwd:/etc/passwd:ro -v /etc/group:/etc/group:ro --gpus=1 --user $(id -u):$(id -g) ontresearch/dorado /bin/bash -c "cd /home/elruiz && dorado basecaller --emit-fastq dna_r10.4.1_e8.2_400bps_sup\@v4.2.0 A_data/merged_pod5/unfiltered_seabass_edna_reads.pod5 > A_data/basecalling_sup/unfiltered_seabass_edna_reads_sup.fastq"





##### PART 3 - CPU server (bash): Installing all required packages for all scripts #####

## Installing micromamba

# sudo apt-get update
# curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
# eval "$(./bin/micromamba shell hook -s posix)"
# ./bin/micromamba shell init -s bash -p ~/micromamba
# micromamba activate
# mamba create -n modbam2bed -c bioconda -c conda-forge -c epi2melabs modbam2bed
# sudo apt-get update


## Installing packages outside with sudo

# sudo apt-get update
# sudo apt-get install python3-dev
# sudo apt-get install samtools
# sudo apt-get install minimap2
# sudo apt-get update


## Installing required packages within a virtual environment

# python3.8 -m venv --prompt remora --copies venv_remora
# source venv_remora/bin/activate
# pip install --upgrade pip setuptools wheel
# pip install --upgrade pip pysam wheel
# pip install --upgrade pip ont-remora wheel
# pip install ont-remora
# pip install pod5
# deactivate


## Installing guppy

# sudo apt-get update
# sudo apt install wget lsb-release
# export PLATFORM=$(lsb_release -cs)
# wget -O- https://cdn.oxfordnanoportal.com/apt/ont-repo.pub | sudo apt-key add -
# echo "deb http://cdn.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
# sudo apt-get update
# sudo apt install ont-guppy


## Installing dorado and basecalling models

# sudo apt-get update
# wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.3.2-linux-x64.tar.gz
# tar -xf /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/dorado-0.3.2-linux-x64.tar.gz
# pod5 inspect debug /mnt/d/Users/'Eliot RUIZ'/Documents/Methylation/A_data/merged_pod5/unfiltered_seabass_edna_reads.pod5 | grep sample_rate # Need to use 5 kHz models
# dorado-0.3.2-linux-x64/bin/dorado download --model all
# sudo apt-get update





##### PART 4 - CPU server (bash): Demultiplexing the sup basecalled FASTQ #####

# SECONDS=0
# guppy_barcoder \
# --input_path /mnt/d/Users/'Eliot RUIZ'/Documents/Methylation/A_data/basecalling_sup \
# --save_path /mnt/d/Users/'Eliot RUIZ'/Documents/Methylation/A_data/demultiplexing_sup \
# --barcode_kits SQK-NBD114-24
# duration=$SECONDS
# echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." ## ABOUT 1H





##### PART 5 - Local computer (R): Converting FASTQ files in FASTA files and saving names #####

## Initialisation

# renv::init()
# renv::snapshot()
# renv::status()

# renv::restore() # Line to run to directly install dependencies of the whole project with the right versions
library(ShortRead)
library(Biostrings)
library(tidyverse)

path_h_drive_a = "H:/seabass_edna_methylation_data/A_data/"
# Note: The full data folder is very heavy so it is only stored on an external hard drive (H) and the CPU server


## Saving all sup basecalled reads in a FASTA file

all_reads = c()
all_pod5_names = tibble()

for(i in 1:24){
  
  cat(paste0("Barcode number", i, "\n"))
  
  path_barcode = paste0(path_h_drive_a, "demultiplexing_sup/barcode", ifelse(i < 10, paste0("0", i), i))
  
  if(length(list.files(path_barcode)) > 0){
    
    fastq_files = list.files(path = path_barcode, pattern = "*.fastq", full.names = T)
    
    for(j in 1:length(fastq_files)){
      
      temp_reads = readDNAStringSet(filepath = fastq_files[j], format = "fastq")
      
      pod5_names = names(temp_reads)
      
      names(temp_reads) = paste0("bar", i, "_seq", 1:length(temp_reads), "_wid", width(temp_reads))
      
      all_pod5_names = rbind(all_pod5_names, tibble(BARCODE = i, INITIAL_NAMES = word(pod5_names, 1), 
                                                    NEW_NAMES = names(temp_reads)))
      
      all_reads = c(all_reads, temp_reads)
      
    }
    
  }
  
}

writeFasta(all_reads[[1]], paste0(path_h_drive_a, "merged_sup/all_sup_eDNA_reads_full.fasta"), mode = "w")
lapply(all_reads[2:length(all_reads)], function(x) 
  writeFasta(x, paste0(path_h_drive_a, "merged_sup/all_sup_eDNA_reads_full.fasta"), mode = "a"))

all_reads_full = readDNAStringSet(paste0(path_h_drive_a, "merged_sup/all_sup_eDNA_reads_full.fasta"), format = "fasta")
all_reads_full

table(word(names(all_reads_full), 1, sep = fixed("_")))


## Saving sup pod5 names along with the new names

write.csv(all_pod5_names, paste0(path_h_drive_a, "merged_sup/all_sup_pod5_names.csv"), row.names = F)
all_pod5_names = tibble(read.csv(paste0(path_h_drive_a, "merged_sup/all_sup_pod5_names.csv")))
all_pod5_names
