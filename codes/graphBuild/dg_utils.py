import numpy as np
import json
import os
import shutil
import random


def ClearDir(dirpath):
    if os.path.exists(dirpath):
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


def GetPointSets(pts_size,
                 num_clusters,
                 half_space_dist_,
                 dim,
                 gauss_delta_,
                 _means=np.array([])):  # 获得最初的point set（包含num_clusters个簇，每个簇
    _points = []
    _labels = []
    # 生成一个多维高斯分布
    for i in range(num_clusters):
        if _means.shape[0] == 0:
            _mean = np.random.uniform(-half_space_dist_, half_space_dist_,
                                      (dim))  # z
        else:
            _mean = _means[i]

        _cov = np.diag(np.array([gauss_delta_ for i in range(dim)]))  # 方差0.5
        _points.append(
            np.random.multivariate_normal(size=pts_size // num_clusters,
                                          mean=_mean,
                                          cov=_cov))
        _labels.append([i for t in range(pts_size // num_clusters)])

    p = np.array(_points).reshape((pts_size, dim))
    l = np.array(_labels).reshape((pts_size, 1))

    state = np.random.get_state()
    np.random.shuffle(p)
    np.random.set_state(state)
    np.random.shuffle(l)
    return p, l


def GenDistubIds(pts_size, keep_ratio):
    ids = range(0, pts_size)

    # keep_ids 相似性的点
    # dist_ids 不相似的点
    keep_ids = random.sample(ids, int(keep_ratio * pts_size))
    keep_ids.sort()
    dist_ids = [i for i in ids if i not in keep_ids]

    return keep_ids, dist_ids

    # keep_ids, dist_ids = GenDistubIds(pts_size, keep_ratio)


# disturb given points
def DisturbPoints(inputs, dim, keep_ids, dist_ids, disturb_dist, HARD_MOVE):
    output = inputs

    for i in dist_ids:
        if HARD_MOVE == True:
            _moveVec = np.tile(disturb_dist, dim)
        else:
            _moveVec = np.random.uniform(-disturb_dist, disturb_dist, (dim))

        output[i] += _moveVec
    return output, keep_ids


# disturb clusters of given label
def DisturbClusters(inputs, dim, labels, disturb_label, pts_size, disturb_dist,
                    HARD_MOVE):
    ids = range(0, pts_size)
    # keep_ids 相似性的点
    # dist_ids 不相似的点
    dist_ids = [i for i in range(len(labels)) if labels[i] == disturb_label]
    keep_ids = [i for i in ids if i not in dist_ids]

    return DisturbPoints(inputs, dim, keep_ids, dist_ids, disturb_dist,
                         HARD_MOVE)


# translate the clusters
def TranslateClusters(inputs, dim, labels, disturb_labels, pts_size,
                      disturb_dist, HARD_MOVE):
    outputs = inputs

    ids = range(0, pts_size)
    # dist_ids 不相似的点
    dist_ids = [i for i in range(len(labels)) if labels[i] in disturb_labels]
    # keep_ids 相似性的点
    keep_ids = [i for i in ids if i not in dist_ids]

    if HARD_MOVE == True:
        _moveVec = np.tile(disturb_dist, dim)
    else:
        _moveVec = np.random.uniform(-disturb_dist, disturb_dist, (dim))

    # the same displacement
    for i in dist_ids:
        outputs[i] += _moveVec
    return outputs, keep_ids


# disturb point sets from different clusters
def DisturbPointSets(inputs, labels, disturb_label_num, keep_ratio):
    assert (disturb_label_num > 0 and disturb_label_num <= num_clusters)

    output = inputs
    ids = range(0, pts_size)
    # keep_ids 相似性的点
    # dist_ids 不相似的点
    disturb_labels = random.sample(range(num_clusters), disturb_label_num)

    num_each_cluster = int((1. - keep_ratio) * pts_size / disturb_label_num)

    # for each cluster disturb the same number of
    dist_ids = {}
    for label in disturb_labels:
        label_ids = [i for i in range(len(labels)) if labels[i] == label]
        dist_ids[label] = random.sample(label_ids, num_each_cluster)

    keep_ids = []
    for i in ids:
        flag = True
        for label in disturb_labels:
            if i in dist_ids[label]:
                flag = False
        if flag == True:
            keep_ids.append(i)

    # ...
    for label in disturb_labels:
        if HARD_MOVE_ == True:
            _moveVec = np.tile(disturb_dist, dim) * random.sample([-1, 1],
                                                                  1)[0]
        else:
            _moveVec = np.random.uniform(-disturb_dist, disturb_dist, (dim))

        # the same displacement
        for i in dist_ids[label]:
            output[i] += _moveVec

    return output, keep_ids


def shiftAllPoints(inputs, dim, disturb_dist, HARD_MOVE_):
    if HARD_MOVE_ == True:
        _moveVec = np.tile(disturb_dist, dim) * random.sample([-1, 1], 1)[0]
    else:
        _moveVec = np.random.uniform(-disturb_dist, disturb_dist, (dim))

    output = inputs
    for i in range(output.shape[0]):
        output[i] += _moveVec

    # no keeping id
    return output, []


# add noise to undisturbed points
def addNoise(inputs, dim, dist_ids, noise_intensity):
    outputs = inputs
    for id in dist_ids:
        moveVec = np.random.uniform(-noise_intensity, noise_intensity, (dim))
        outputs[id] += moveVec

    print("add noise to current data")
    return outputs


def overlapClusters(inputs, pts_size, dim, means, labels, merge_labels):
    output = inputs
    ids = range(0, pts_size)

    merge_mean = np.zeros((dim))
    # 计算出要合并的几个簇中心的中心
    for label in merge_labels:
        merge_mean += means[label]
    merge_mean /= len(merge_labels)

    for id in ids:
        # move the cluster center to the same location
        if labels[id] in merge_labels:
            # print(output[id, :].shape)
            # print(merge_mean.shape)
            # print(means[labels[id], :].shape)
            output[id, :] += merge_mean - means[labels[id], :].reshape(dim, )

    return output


# 缩小整个类簇而不改变knn
def scaleCluster(inputs,
                 pts_size,
                 dim,
                 labels,
                 shrink_label,
                 new_center,
                 scale_factor=0.25):
    output = inputs
    ids = range(0, pts_size)

    # 先算出该类簇中心
    scale_center = cluster_center(inputs=inputs,
                                  labels=labels,
                                  label=shrink_label,
                                  pts_size=pts_size,
                                  dim=dim)

    # 将所有点围绕该中心进行缩放,并移动到
    shrink_ids = [id for id in ids if labels[id] == shrink_label]
    for id in shrink_ids:
        output[id] = new_center + (output[id] - scale_center) * scale_factor

    return output


def cluster_center(inputs, pts_size, dim, labels, label):
    ids = range(0, pts_size)

    center = np.zeros((dim))
    count = 0
    for id in ids:
        if labels[id] == label:
            center += inputs[id, :]
            count += 1

    center /= float(count)
    return center


def splitClusters(inputs, pts_size, dim, labels, split_label, disturb_dist):
    output = inputs
    ids = range(0, pts_size)

    split_ids = [id for id in ids if labels[id] == split_label]

    split_ids_0 = random.sample(split_ids, int(len(split_ids) / 2))
    split_ids_1 = [id for id in split_ids if id not in split_ids_0]

    _moveVec = np.tile(-disturb_dist, dim)
    for id in split_ids_0:
        output[id] += _moveVec

    _moveVec = np.tile(disturb_dist, dim)
    for id in split_ids_1:
        output[id] += _moveVec

    return output


def DistOfEdges(dists, indices):
    E = {}
    for i in range(indices.shape[0]):
        vi = i
        for j in range(indices[i].shape[0]):
            vj = indices[i][j]
            if vi != vj:
                E[(vi, vj)] = dists[i][j]
    return E


# save points and labels into file
def savetxt(filepath, cur_points, labels):
    with open(filepath, 'w') as f:
        for i in range(cur_points.shape[0]):
            fstr = ""
            for j in range(cur_points.shape[1]):
                fstr += "%.16f\t" % (cur_points[i][j])
            fstr += "%d\n" % (labels[i])
            f.write(fstr)

    print(os.path.abspath(filepath) + " saved.")


def queryLabelIds(ids, labels, queryLabels=[]):
    return [id for id in ids if labels[id] in queryLabels]


def DistOfEdges(dists, indices):
    E = {}
    for i in range(indices.shape[0]):
        vi = i
        for j in range(indices[i].shape[0]):
            vj = indices[i][j]
            if vi != vj:
                E[(vi, vj)] = dists[i][j]
    return E


def writeInfo(filepath, info):
    with open(filepath, 'w', encoding='utf-8') as f:
        json.dump(info, f)


# # Amplify the cluster
# def DiffuseCluster(inputs, labels, disturb_label):
#     output = inputs
#     assert(disturb_label>=0 and disturb_label<=num_clusters)

#     # compute the centroid of the cluster
#     disturb_ids = [i for i in range(len(labels)) if labels[i] == disturb_label]
#     disturb_inputs = [inputs[id] for id in disturb_ids]
#     disturb_inputs = np.array(disturb_ids)

#     centroid = np.sum(disturb_ids, 0) / disturb_ids.shape(0)

#     # for each point, compute the vector
#     for i in disturb_inputs.shape[0]:
#         # each point go through the vector
#         output[i] +=

#     return output

# # random select keeping edges. NEVER USE IT
# def DisturbEdges(inputs, kd_tree, keep_ratio=0.7):
#     outputs = inputs
#     dists, indices = kd_tree.query(
#         inputs, k=k_closest_count)  # 一口气对所有points构建knn
#     edge_size = indices.shape[0]*(k_closest_count-1)

#     np.arange(0, points.shape[0])
#     ids = range(0, pts_size)

#     # keep_edges 相似性的边
#     keep_edges = []
#     while len(keep_edges) != edge_size*keep_ratio:
#         keep_ids_0 = np.random.randint(0, pts_size)
#         keep_ids_1 = np.random.randint(0, pts_size)
#         # ensure no self-loop and duplicate
#         if keep_ids_0 == keep_ids_1 or [keep_ids_0, keep_ids_1] in keep_edges:
#             continue
#         keep_edges.append([keep_ids_0, keep_ids_1])

#     keep_edges.sort()
#     # dist_edges 不相似的边
#     dist_edges = []
#     for i in range(len(indices)):
#         for j in indices[i]:
#             # print(i)
#             # print(j)
#             if [i, j] not in keep_edges:
#                 dist_edges.append([i, j])

#     # disturb edge endpoints
#     for [i, j] in dist_edges:
#         _moveVec_i = np.random.uniform(-0.5, 0.5, (dim))
#         _moveVec_j = np.random.uniform(-0.5, 0.5, (dim))

#         outputs[i] += _moveVec_i
#         outputs[j] += _moveVec_j

#     return outputs, keep_edges
