#!/bin/bash
# remove generated python files

find . -type d -name '__pycache__' -prune -exec rm -rf {} \;
