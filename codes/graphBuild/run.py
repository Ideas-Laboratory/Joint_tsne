# encoding = utf-8
import os
import sys
import numpy as np
from sklearn.neighbors import KDTree
import json
import time

dirname = os.path.dirname(__file__)

sys.path.append(os.path.join(dirname, ".."))
from thesne.sne_config import SNEConfig
sys.path.append(os.path.join(dirname, "../dataGen/"))
import dg_utils


# read high dimensional dataset from file
def ReadPointSets(filepath, norm_data):
    _points = []
    _labels = []

    _dim = 0
    _pts_size = 0
    _flag = False

    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        _pts_size = len(lines)

        for line in lines:
            items = line.split()
            _point = []

            for item in items[:-1]:
                _point.append(float(item))

                if not _flag:
                    _dim += 1
            _flag = True

            _points.append(_point)
            _labels.append(items[-1])

    p = np.array(_points)
    if norm_data == True:
        min_p = np.min(p)
        max_p = np.max(p)
        p = (p - min_p) / (max_p - min_p)

    return p, np.array(_labels), _pts_size, _dim


def sav_graph(filepath, points, labels, k_closest_count):
    tree = KDTree(points)
    dists, indices = tree.query(points,
                                k=k_closest_count + 1)  # build knn graph from points

    with open(filepath, 'w') as file:
        # write in the number of nodes
        file.write(str(points.shape[0]) + "\n")

        # write each point
        for i in range(points.shape[0]):
            file.write(str(i) + "\t")
            count = 0
            for t in range(indices[i].shape[0]):
                if indices[i][t] == i:
                    continue
                if count == k_closest_count:  # write the k_closest_count neighbors
                    break
                # write incient point and corresponding distance
                file.write(str(indices[i][t]) + "\t" + str(dists[i][t]) + "\t")
                count += 1
            file.write('\n')


if __name__ == '__main__':
    argv = sys.argv
    assert (len(argv) == 2)
    config_path = argv[1]

    start = time.perf_counter()
    with open(config_path, 'r') as f:
        config_json = json.loads(f.read())
        config = SNEConfig(config_json)

        input_dir = config.input_dir
        graph_path = config.graph_dir
        dg_utils.ClearDir(graph_path)

        k_closest_count = config.k_closest_count

        raw_files = []
        for filename in os.listdir(input_dir):
            raw_file = os.path.join(input_dir, filename)
            (file, ext) = os.path.splitext(raw_file)
            if ext == ".txt":
                raw_files.append(raw_file)

        for filepath in raw_files:
            data_id = dg_utils.GetGraphIDFromPath(filepath)
            print("Processing graph " + str(data_id))

            cur_points, labels, pts_size, dim = ReadPointSets(filepath,
                                                              norm_data=True)
            print((pts_size, dim))

            sav_graph(os.path.join(graph_path, "g_{}.txt".format(data_id)),
                      cur_points, labels, k_closest_count)

    elapsed = (time.perf_counter() - start)
    print("Total time for building knn graph:", elapsed)