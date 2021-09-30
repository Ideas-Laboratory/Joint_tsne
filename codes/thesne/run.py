#!/usr/bin/env python3
# encoding=utf-8
import os
import json
import numpy as np

import time
import tsne
import joint_tsne

import sys
sys.path.append("..")
from thesne.sne_config import SNEConfig
import utils.dr_utils as dr_utils

if __name__ == "__main__":
    argv = sys.argv
    assert (len(argv) == 2)
    config_path = argv[1]

    config_json = json.loads(open(config_path).read())
    config = SNEConfig(config_json)

    # For all the methods, the first frame has the same initial embeddings
    dr_utils.ClearDir(config.output_root_dir)

    # initial layout
    y_init = np.random.randn(config.pts_size, 2)
    joint_tsne.run(config, y_init)
