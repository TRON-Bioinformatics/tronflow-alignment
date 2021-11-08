#!/bin/bash

output_folder=output/test3
nextflow main.nf -profile test,conda --library single --output $output_folder
test -s $output_folder/TESTX_S1_L001.bam || { echo "Missing test 3 output file!"; exit 1; }
test -s $output_folder/TESTX_S1_L002.bam || { echo "Missing test 3 output file!"; exit 1; }