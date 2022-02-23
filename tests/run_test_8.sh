#!/bin/bash

output_folder=output/test8
nextflow main.nf -profile test,conda --algorithm mem --skip_trimming --output $output_folder
test -s $output_folder/TESTX_S1_L001.bam || { echo "Missing test 8 output file!"; exit 1; }
test -s $output_folder/TESTX_S1_L001.bam.bai || { echo "Missing test 8 output file!"; exit 1; }
if [ -f $output_folder/TESTX_S1_L001.fastp_stats.html ]; then
    echo "This file should not exist!"; exit 1;
fi
if [ -f $output_folder/TESTX_S1_L001.fastp_stats.json ]; then
    echo "This file should not exist!"; exit 1;
fi

test -s $output_folder/TESTX_S1_L002.bam || { echo "Missing test 8 output file!"; exit 1; }
test -s $output_folder/TESTX_S1_L002.bam.bai || { echo "Missing test 8 output file!"; exit 1; }
if [ -f $output_folder/TESTX_S1_L002.fastp_stats.html ]; then
    echo "This file should not exist!"; exit 1;
fi
if [ -f $output_folder/TESTX_S1_L002.fastp_stats.json ]; then
    echo "This file should not exist!"; exit 1;
fi