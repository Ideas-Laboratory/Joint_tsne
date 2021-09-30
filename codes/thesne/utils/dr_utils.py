import matplotlib.pyplot as plt
import numpy as np
import os
import shutil

def create_juxtaposed_figure(Y_1,
                             Y_2,
                             labels,
                             title_1,
                             title_2,
                             cmap,
                             sub_indices=[],
                             fig_path=None):
    plt.figure(figsize=(12, 4))

    plt.subplot(1, 2, 1)
    plt.title(title_1)
    plt.scatter(Y_1[:, 0], Y_1[:, 1], c=labels, cmap=cmap)

    plt.subplot(1, 2, 2)
    plt.title(title_2)
    plt.scatter(Y_2[:, 0], Y_2[:, 1], c=labels, cmap=cmap)

    # highlight common points
    if len(sub_indices) > 0:
        for i in range(Y_2.shape[0]):
            if i not in sub_indices:
                plt.scatter(Y_2[i, 0], Y_2[i, 1], c='r')

    #
    if fig_path is not None:
        plt.savefig(fig_path, bbox_inches='tight')


def trivial_point_sims(N):
    return {i: 1.0 for i in range(N)}


def select_anchor_point(match_points, rate, pick_from_high=True):
    sim_scores = list(match_points.values())

    if pick_from_high == True:
        # 从前往后是similarity score减小的下标
        rank = np.argsort(sim_scores)[::-1]
    else:
        rank = np.argsort(sim_scores)

    points = list(match_points.keys())
    anchor_point = {
        points[i]: sim_scores[i]
        for i in rank[:int(rate * len(rank))]
    }

    return anchor_point


# select edge with lowest similarity
def select_anchor_edge(match_edges, rate, pick_from_high=True):
    sim_scores = list(match_edges.values())

    if pick_from_high == True:
        # 从前往后是similarity score变大的下标
        rank = np.argsort(sim_scores)[::-1]
    else:
        rank = np.argsort(sim_scores)

    edges = list(match_edges.keys())
    anchor_edges = {
        edges[i]: sim_scores[i]
        for i in rank[:int(rate * len(rank))]
    }

    return anchor_edges


# compute edge similarity based on point similarity
def find_match_edges(edges0, edges1, match_points):
    match_edges = {}
    for e0 in edges0:
        for e1 in edges1:
            if e0 == e1:  # common edge
                match_edges[e0] = match_points[e0[0]] * match_points[e0[1]]
                # match_edges[(e0, e1)] = 1 - match_points[(e0[0], e0[0])]*match_points[(e0[1], e0[1])]
                #(match_points[(e0[0], e0[0])] + match_points[(e0[1], e0[1])] )/2.

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
    # print("read file: " + filepath)

    match_points = {}
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip('\n')
            items = line.split()
            match_points[int(items[0])] = float(items[1])

    return match_points


# preload synset, only return the serial number
def load_synset(file_path="../dataGen/vgg16-tf-machrisaa/synset.txt"):
    synset = {}
    for l in open(file_path).readlines():
        ss = l.strip().split(maxsplit=1)
        code = ss[0]
        name = ss[1].strip().split(',')[0]
        synset[code] = name

    print("synset loaded.")
    return synset


def get_label_set(labels):
    labelSet = []
    for label in labels:
        if label not in labelSet:
            labelSet.append(label)

    labelSet.sort()
    return labelSet


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
        # print("你真的看清楚你要干什么了吗?")
        # input()
        print("正在删除.....", dirpath)
        shutil.rmtree(path=dirpath)
    os.makedirs(dirpath)


def GetFilesIn(dir):
    res = []
    for dir, _, files in os.walk(dir):
        for file in files:
            res.append(dir + file)
    return res


def GetGraphIDFromPath(path):
    print(path)
    ID = path.split(".")[-2].split("_")[-1]
    return int(ID)


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