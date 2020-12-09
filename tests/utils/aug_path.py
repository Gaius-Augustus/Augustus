#!/usr/bin/env python3

import os
import shutil

def rmtree_if_exists(dir, even_none_empty=False):
    if os.path.exists(dir):
        if even_none_empty or len(os.listdir(dir)) == 0:
            shutil.rmtree(dir)


def mkdir_if_not_exists(dir):
    if not os.path.exists(dir):
        os.makedirs(dir, exist_ok=True)
