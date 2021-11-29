import matplotlib.pyplot as plt
import numpy as np
import os
import shutil

# compute edge similarity based on point similarity
def find_match_edges(edges0, edges1, match_points):
    match_edges = {}
    for e0 in edges0:
        for e1 in edges1:
            if e0 == e1:  # common edge
                match_edges[e0] = match_points[e0[0]] * match_points[e0[1]]
    return match_edges


def common_neighbor_ratio(node_num, edges0, edges1, Neighbor_size):
    ratios = [0 for i in range(node_num)]
    for e0 in edges0:
        for e1 in edges1:
            if e0 == e1:
                ratios[e0[0]] += 1

    for i in range(node_num):
        ratios[i] /= float(Neighbor_size)

    return ratios


def read_match_points(filepath):
    match_points = {}
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip('\n')
            items = line.split()
            match_points[int(items[0])] = float(items[1])

    return match_points


# save points and labels into file
def saveVectors(filepath, cur_points, labels):
    with open(filepath, 'w') as f:
        for i in range(cur_points.shape[0]):
            fstr = ""
            for j in range(cur_points.shape[1]):
                fstr += "%.16f\t" % (cur_points[i][j])
            fstr += "%s\n" % (labels[i])
            # fstr += str(labels[i])
            f.write(fstr)

    print(filepath + " saved.")


def ClearDir(dirpath):
    if os.path.exists(dirpath):
        print("Deleting", dirpath)
        shutil.rmtree(path=dirpath)
    os.makedirs(dirpath)


def read_graph(filepath, K):
    edges = []
    # here we continue to readline
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            line = line.strip('\n')
            items = line.split()

            # discard the first line
            if len(items) == 1:
                continue

            # print(items)
            assert (len(items) == K * 2 + 1)

            i = int(items[0])
            # WE discard distance here
            for n in range(K):
                j = int(items[n * 2 + 1])
                # n*2+2
                edges.append((i, j))

    return edges


def norm_data(X):
    min_x = np.min(X)
    max_x = np.max(X)
    X = (X - min_x) / (max_x - min_x)
    print("Normalize data with: {}, {}".format(min_x, max_x))
    return X