#!/bin/bash

output_folder=output/test7
nextflow main.nf -profile test,conda --output $output_folder --input_files false \
--input_fastq1 test_data/TESTX_S1_L001_R1_001.fastq.gz \
--library single --input_name test
test -s $output_folder/test.bam || { echo "Missing test 7 output file!"; exit 1; }