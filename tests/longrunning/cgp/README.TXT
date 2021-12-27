# Test of comparative gene prediction on large realistic dataset that allows a monitoring of the accuracy.

### Requirements:

1) Install required python packages
pip3 install wget gitpython psutil memory-profiler

2) eval (Evan Keibler)
Copy perl sources to a dir $EVAL_DIR of your choosing, e.g.
EVAL_DIR=/home/mario/tools/eval-2.2.8


### Run with

./executeTestCGP.py --chunks 27 30 47 54 57 80 86 101 118 --eval --augustusDir ../../../ --dataDir data --evalDir $EVAL_DIR --jobs 5

Further options:
--pathToGitRepo: Path to the Augustus Git repository. If the
        path is specified, the file additional_information.json
        is supplemented with information from the Git repository,
        such as the current revisoin hash or the commit date.
        (Not recommended for local execution.)

Output:
The evaluation output can be found in the folder
cgp/data/cgp12way/ACCURACY.
Furthermore, the file additional_information.json contains
information about the runtime or the memory usage as well as
git repository information if the option --pathToGitRepo is set.

### Extract desired accuracy values and save them as JSON file
(execution from the directory tests/longrunning/)

./extract_accuracy_values.py --path cgp/data/cgp12way/ACCURACY
