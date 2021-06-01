
all : clean test check

clean:
	rm -rf output
	rm -f .nextflow.log*
	rm -rf .nextflow*

test:
	nextflow main.nf --help
	nextflow main.nf -profile test,conda --output output/test1
	nextflow main.nf -profile test,conda --inception --output output/test2
	nextflow main.nf -profile test,conda --library single --output output/test3
	nextflow main.nf -profile test,conda --algorithm mem --output output/test4
	nextflow main.nf -profile test,conda --algorithm mem --library single --output output/test5
	nextflow main.nf -profile test,conda --output output/test6 --input_files false \
	--input_fastq1 test_data/TESTX_S1_L001_R1_001.fastq.gz \
	--input_fastq2 test_data/TESTX_S1_L001_R2_001.fastq.gz --input_name test
	nextflow main.nf -profile test,conda --output output/test7 --input_files false \
	--input_fastq1 test_data/TESTX_S1_L001_R1_001.fastq.gz \
	--library single --input_name test

check:
	test -s output/test1/TESTX_S1_L001.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test1/TESTX_S1_L002.bam || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test2/TESTX_S1_L001.bam || { echo "Missing test 2 output file!"; exit 1; }
	test -s output/test2/TESTX_S1_L002.bam || { echo "Missing test 2 output file!"; exit 1; }
	test -s output/test3/TESTX_S1_L001.bam || { echo "Missing test 3 output file!"; exit 1; }
	test -s output/test3/TESTX_S1_L002.bam || { echo "Missing test 3 output file!"; exit 1; }
	test -s output/test4/TESTX_S1_L001.bam || { echo "Missing test 4 output file!"; exit 1; }
	test -s output/test4/TESTX_S1_L002.bam || { echo "Missing test 4 output file!"; exit 1; }
	test -s output/test5/TESTX_S1_L001.bam || { echo "Missing test 5 output file!"; exit 1; }
	test -s output/test5/TESTX_S1_L002.bam || { echo "Missing test 5 output file!"; exit 1; }
	test -s output/test6/test.bam || { echo "Missing test 6 output file!"; exit 1; }
	test -s output/test7/test.bam || { echo "Missing test 7 output file!"; exit 1; }
