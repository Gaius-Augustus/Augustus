#### Requires:

1) GNU Parallel. Install with
sudo apt install parallel

2) Python Memory Profiler. Install with
pip install memory-profiler

3) eval (Evan Keibler)
Copy perl sources to a dir $EVAL_DIR of your choosing, e.g.
EVAL_DIR=/home/mario/tools/eval-2.2.8

### Run with

./test_single.py --pathToGitRepo ../../../.git --evalDir $EVAL_DIR --jobs 8
