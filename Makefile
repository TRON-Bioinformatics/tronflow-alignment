clean:
	rm -rf output
	rm -rf work
	rm -f report.html*
	rm -f timeline.html*
	rm -f trace.txt*
	rm -f dag.dot*
	rm -f .nextflow.log*
	rm -rf .nextflow*

build-docker:
	docker build -t tron-bioinformatics/tronflow-bwa:1.1.0 .

test-docker:
	nextflow main.nf --input_files test_data/test_input_paired.txt -profile docker --reference `pwd`/test_data/ucsc.hg19.minimal.fasta

test-conda-aln-paired:
	nextflow main.nf --input_files test_data/test_input_paired.txt -profile conda --reference `pwd`/test_data/ucsc.hg19.minimal.fasta

test-conda-aln-single:
	nextflow main.nf --input_files test_data/test_input_single.txt -profile conda --reference `pwd`/test_data/ucsc.hg19.minimal.fasta --library single

test-conda-mem-paired:
	nextflow main.nf --input_files test_data/test_input_paired.txt -profile conda --reference `pwd`/test_data/ucsc.hg19.minimal.fasta --algorithm mem

test-conda-mem-single:
	nextflow main.nf --input_files test_data/test_input_paired.txt -profile conda --reference `pwd`/test_data/ucsc.hg19.minimal.fasta --algorithm mem --library single
