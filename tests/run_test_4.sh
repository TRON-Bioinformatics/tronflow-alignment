#!/bin/bash

output_folder=output/test4
nextflow main.nf -profile test,conda --algorithm mem --output $output_folder
test -s $output_folder/TESTX_S1_L001.bam || { echo "Missing test 4 output file!"; exit 1; }
test -s $output_folder/TESTX_S1_L002.bam || { echo "Missing test 4 output file!"; exit 1; }