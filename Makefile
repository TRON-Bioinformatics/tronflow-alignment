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
	docker build -t tron-bioinformatics/tronflow-bwa:1.0.0 .

test-docker:
	/home/priesgo/bin/nextflow main.nf --input_files test_data/test_input.txt -profile docker --reference `pwd`/test_data/ucsc.hg19.minimal.fasta

test-conda:
	nextflow main.nf --input_files test_data/test_input.txt -profile conda --reference `pwd`/test_data/ucsc.hg19.minimal.fasta
