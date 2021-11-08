#!/bin/bash

output_folder=output/test1
nextflow main.nf -profile test,conda --output $output_folder
test -s $output_folder/TESTX_S1_L001.bam || { echo "Missing test 1 output file!"; exit 1; }
test -s $output_folder/TESTX_S1_L002.bam || { echo "Missing test 1 output file!"; exit 1; }