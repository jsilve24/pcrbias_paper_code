#!/bin/bash
#SBATCH -J demultiplex
#SBATCH --mem=64000
#SBATCH --output=demultiplex.out
#SBATCH --error=demultiplex.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=justin.silverman@duke.edu

PATH1=2_sync_barcodes/run1
MAP1=0_mapping/run1/2017.08.03MappingFile_corrected.txt

module load qiime

# Demultiplexing for DADA2 import that Qiime does not do any quality filtering or else forward and reverse read order may not match exactly
# http://benjjneb.github.io/dada2/faq.html

# Split_Libraries_Fastq
# First split libraries (just modifies reads/header by matching sampleID from mapping file to barcodes) note some reads are lost because of barcode error. 
split_libraries_fastq.py --store_demultiplexed_fastq --rev_comp_mapping_barcodes -r 999 -n 999 -q 0 -p 0.0001 -i $PATH1/R1.paired.fastq.gz -b $PATH1/I1.synced.fastq.gz -m $MAP1 -o 3_demultiplex/s1r1_split_libraries/
split_libraries_fastq.py --store_demultiplexed_fastq --rev_comp_mapping_barcodes -r 999 -n 999 -q 0 -p 0.0001 -i $PATH1/R2.paired.fastq.gz -b $PATH1/I1.synced.fastq.gz -m $MAP1 -o 3_demultiplex/s1r2_split_libraries/


# Split_Sequence_File_On_Sample_Ids
# Next take each of those seqs.fastq files and break them apart into distinct files based on the header information
split_sequence_file_on_sample_ids.py -i 3_demultiplex/s1r1_split_libraries/seqs.fastq --file_type fastq -o 3_demultiplex/s1r1_split_samples/
split_sequence_file_on_sample_ids.py -i 3_demultiplex/s1r2_split_libraries/seqs.fastq --file_type fastq -o 3_demultiplex/s1r2_split_samples/

# Gzip results
gzip 3_demultiplex/s1r1_split_samples/*.fastq
gzip 3_demultiplex/s1r2_split_samples/*.fastq


