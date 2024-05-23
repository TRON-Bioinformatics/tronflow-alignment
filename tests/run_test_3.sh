#!/bin/bash

output_folder=output/test3
nextflow main.nf -profile test,conda,ci --library single --output $output_folder
test -s $output_folder/TESTX_S1_L001.bam || { echo "Missing test 3 output file!"; exit 1; }
test -s $output_folder/TESTX_S1_L001.bam.bai || { echo "Missing test 3 output file!"; exit 1; }
test -s $output_folder/TESTX_S1_L001.fastp_stats.html || { echo "Missing test 3 output file!"; exit 1; }
test -s $output_folder/TESTX_S1_L001.fastp_stats.json || { echo "Missing test 4 output file!"; exit 1; }
test -s $output_folder/TESTX_S1_L002.bam || { echo "Missing test 3 output file!"; exit 1; }
test -s $output_folder/TESTX_S1_L002.bam.bai || { echo "Missing test 3 output file!"; exit 1; }
test -s $output_folder/TESTX_S1_L002.fastp_stats.html || { echo "Missing test 3 output file!"; exit 1; }
test -s $output_folder/TESTX_S1_L002.fastp_stats.json || { echo "Missing test 3 output file!"; exit 1; }