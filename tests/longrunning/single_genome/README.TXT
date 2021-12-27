#### Requires:

1) Install required python packages
pip3 install wget gitpython psutil memory-profiler

2) eval (Evan Keibler)
Copy perl sources to a dir $EVAL_DIR of your choosing, e.g.
EVAL_DIR=/home/mario/tools/eval-2.2.8

### Run with

./test_single.py --evalDir $EVAL_DIR --jobs 9

Further options:
--keepData: If this option is used, the downloaded data will not
        be deleted and it is not necessary to download it again.
--pathToGitRepo: Path to the Augustus Git repository. If the
        path is specified, the file additional_information.json
        is supplemented with information from the Git repository,
        such as the current revisoin hash or the commit date.
        (Not recommended for local execution.)

Output:
The output of each AUGUSTUS run can be found in the output folder.
The corresponding evaluation can be found in the folder output/eval.
Furthermore, the file additional_information.json contains information
about the runtime or the memory usage of each run as well as git
repository information if the option --pathToGitRepo is set.


### Extract desired accuracy values and save them as JSON file
(execution from the directory tests/longrunning/)

./extract_accuracy_values.py --path single_genome/output/eval
