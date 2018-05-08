#! /bin/sh
# Test running utrrnaseq via command line call
# Ingo Bulla
# 19 Jan 2018

oneTimeSetUp() { 
	prog=../Debug/utrrnaseq
	test_data_dir=../for-test-cases 
}


oneTimeTearDown() {
	rm utrs.gff
}


# Test standard scenario
testStandard() {
	dir=${test_data_dir}/standard-test-set
	eval $prog --in-scaffold-file $dir/test.fa \
    	       --in-coding-region-file $dir/test-start-stop.gtf \
    	       --in-intron-file $dir/test-introns.hints \
               --in-repeat-file $dir/test-repeat.hints \
               --in-wiggle-file $dir/test.wig
	assertEquals $? 0
}

# Test softmasking scenario, i.e. no repeat file is provided
testSoftmasing() {
	dir=${test_data_dir}/softmasking-test-set
	eval $prog --in-scaffold-file $dir/test.fa \
    	       --in-coding-region-file $dir/test-start-stop.gtf \
    	       --in-intron-file $dir/test-introns.hints \
               --in-wiggle-file $dir/test.wig
	assertEquals $? 0
}

# Test invalid input, i.e. a filename passed by command line does not refer to an existing file
testInvalidFilename() {
	dir=${test_data_dir}/standard-test-set
	eval $prog --in-scaffold-file $dir/test.fa \
    	       --in-coding-region-file $dir/test-start-stop.gtf \
    	       --in-intron-file $dir/test-introns.hints \
               --in-repeat-file $dir/test-repeat.hints \
               --in-wiggle-file $dir/file_does_not_exist.wig
	assertNotEquals $? 0
}


# Load shUnit2.
. shunit2




