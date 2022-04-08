
all : clean test

clean:
	rm -rf output
	rm -f .nextflow.log*
	rm -rf .nextflow*

test:
	bash tests/run_test_0.sh
	bash tests/run_test_1.sh
	bash tests/run_test_2.sh
	bash tests/run_test_3.sh
	bash tests/run_test_4.sh
	bash tests/run_test_5.sh
	bash tests/run_test_6.sh
	bash tests/run_test_7.sh
	bash tests/run_test_8.sh
	bash tests/run_test_9.sh
	bash tests/run_test_10.sh
	# STAR indices take over 1 GB and we did not manage to make it work in GitHub actions
	#bash tests/run_test_11.sh
	#bash tests/run_test_12.sh
