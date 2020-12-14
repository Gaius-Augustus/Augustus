# Test of comparative gene prediction on large realistic dataset that allows a monitoring of the accuracy.

### Requirements:

1) Install required python packages
pip3 install wget gitpython psutil memory-profiler

2) eval (Evan Keibler)
Copy perl sources to a dir $EVAL_DIR of your choosing, e.g.
EVAL_DIR=/home/mario/tools/eval-2.2.8


### Run with

./executeTestCGP.py --chunks 27 30 47 54 57 80 86 101 118 --eval --augustusDir ../../../ --dataDir data --evalDir $EVAL_DIR --jobs 5


### Extract desired accuracy values and save them as JSON file
(execution from the directory tests/longrunning/)

./extract_accuracy_values.py --path cgp/data/cgp12way/ACCURACY
