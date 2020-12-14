#### Requires:

1) GNU Parallel. Install with
sudo apt install parallel

2) Python Memory Profiler. Install with
pip3 install memory-profiler wget gitpython

3) eval (Evan Keibler)
Copy perl sources to a dir $EVAL_DIR of your choosing, e.g.
EVAL_DIR=/home/mario/tools/eval-2.2.8

### Run with

./test_single.py --evalDir $EVAL_DIR --jobs 9

### Extract desired accuracy values and save them as JSON file
(execution from the directory tests/longrunning/)

./extract_accuracy_values.py --path single_genome/output/eval
