#!/bin/bash

output_folder=output/test6
nextflow main.nf -profile test,conda --output $output_folder --input_files false \
--input_fastq1 test_data/TESTX_S1_L001_R1_001.fastq.gz \
--input_fastq2 test_data/TESTX_S1_L001_R2_001.fastq.gz --input_name test
test -s $output_folder/test.bam || { echo "Missing test 6 output file!"; exit 1; }