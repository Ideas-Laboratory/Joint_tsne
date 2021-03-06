import os

import json
import numpy as np
import time

from functools import reduce

from tsne import x2p, pca, tsne

import sys
sys.path.append("..")
import utils.dr_utils as dr_utils
from thesne.sne_config import SNEConfig


def joint_tsne(gamma_,
               Y_0=np.array([]),
               X_1=np.array([]),
               Y_1_I=np.array([]),
               match_edges={},
               no_dims=2,
               initial_dims_1=50,
               perplexity=30.0,
               verbose=1):
    """
        Runs t-SNE on two datasets in the NxD array X to reduce its
        dimensionality to no_dims dimensions. The syntaxis of the function is
        `Y = tsne.tsne(X, no_dims, perplexity), where X is an NxD NumPy array.
    """

    # Check inputs
    if isinstance(no_dims, float):
        print("Error: array X should have type float.")
        return -1
    if round(no_dims) != no_dims:
        print("Error: number of dimensions should be an integer.")
        return -1

    # Initialize variables
    X_1 = pca(X_1, initial_dims_1).real

    (n1, d1) = X_1.shape

    max_iter = 2000
    initial_momentum = 0.5
    final_momentum = 0.8
    # momentum
    eta = 500
    min_gain = 0.01

    Y1 = Y_1_I
    # dy1: gradient for data2
    dY1 = np.zeros((n1, no_dims))
    iY1 = np.zeros((n1, no_dims))
    gains1 = np.ones((n1, no_dims))

    # 对称化
    P1 = x2p(X_1, 1e-5, perplexity)
    P1 = P1 + np.transpose(P1)
    P1 = P1 / np.sum(P1)

    early_exagger = True
    # early exaggeration
    if early_exagger:
        P1 = P1 * 4.
    P1 = np.maximum(P1, 1e-12)

    # acclerate computation
    match_edges_len = len(match_edges)

    # Run iterations
    for iter in range(max_iter):
        # Compute pairwise affinities
        sum_Y1 = np.sum(np.square(Y1), 1)
        num1 = -2. * np.dot(Y1, Y1.T)
        num1 = 1. / (1. + np.add(np.add(num1, sum_Y1).T, sum_Y1))  # qij分子 ????
        num1[range(n1), range(n1)] = 0.
        Q1 = num1 / np.sum(num1)  # qij
        Q1 = np.maximum(Q1, 1e-12)
        ''' 
        Compute gradient
        '''
        PQ1 = P1 - Q1

        for i in range(n1):
            dY1[i, :] = np.sum(
                np.tile(PQ1[:, i] * num1[:, i],
                        (no_dims, 1)).T * (Y1[i, :] - Y1), 0)

        # FIXME: use array instead of dictionary
        for (vi, vj) in match_edges:
            # get edge vectors
            d0 = np.subtract(Y_0[vi, :], Y_0[vj, :])
            d1 = np.subtract(Y1[vi, :], Y1[vj, :])
            tmp = (gamma_ * match_edges[(vi, vj)] * (d0 - d1) /
                   match_edges_len)
            dY1[vi, :] -= tmp
            dY1[vj, :] += tmp

        # Perform the update
        if iter < 20:
            momentum = initial_momentum
        else:
            momentum = final_momentum

        gains1 = (gains1 + 0.2) * ((dY1 > 0.) != (iY1 > 0.)) + \
            (gains1 * 0.8) * ((dY1 > 0.) == (iY1 > 0.))
        gains1[gains1 < min_gain] = min_gain

        # update momentum
        iY1 = momentum * iY1 - eta * (gains1 * dY1)

        # update y
        Y1 = Y1 + iY1

        # subtract mean
        Y1 = Y1 - np.tile(np.mean(Y1, 0), (n1, 1))
        ''' 
        Compute current value of cost function
        '''
        if (iter + 1) % 100 == 0:
            C0 = np.sum(P1 * np.log(P1 / Q1))
            C1 = reduce(
                lambda x, eij: x + match_edges[eij] * np.sum(
                    np.square(
                        np.subtract(
                            np.subtract(Y_0[eij[0], :], Y_0[eij[1], :]),
                            np.subtract(Y1[eij[0], :], Y1[eij[1], :])))) /
                match_edges_len, match_edges, 0)

            if verbose:
                C = C0 + gamma_ * C1
                print(
                    "Iteration %d: KL error is %f, edge error is %f, total error is %f"
                    % (iter + 1, C0, C1, C))

        # Stop lying about P-values
        if early_exagger:
            if iter == 100:
                P1 = P1 / 4.

    # Return solution
    return Y1


def run(config, Y_I, _Y_0_=np.array([])):
    print("Run joint-tsne.")

    # algorithm parameters
    perplexity = config.perplexity
    gamma = config.gamma 

    data_ids = config.data_ids
    data_dims = config.data_dims

    # for each parameter
    param_name = "perp:{}, gamma:{}".format(perplexity, gamma)
    print("parameters:   " + param_name)

    # use t-sne to project the first frame
    data_id_0 = data_ids[0]
    data_dim_0 = data_dims[0]

    data_path = os.path.join(config.input_dir, "f_{}.txt".format(data_id_0))
    X0 = np.loadtxt(data_path, usecols=range(0, data_dim_0))
    if config.norm:
        X0 = dr_utils.norm_data(X0)
    labels = np.loadtxt(data_path, usecols=(data_dim_0, ), dtype=str)

    # if first frame is provided
    if _Y_0_.shape[0] != 0:
        print("use first frame from dynamic t-sne.")
        Y_0 = _Y_0_
    else:
        print("generate first frame from joint t-sne.")
        Y_0, dump = tsne(X=X0,
                         Y_I=Y_I,
                         no_dims=2,
                         initial_dims=data_dim_0,
                         perplexity=perplexity)
    dr_utils.saveVectors(
        os.path.join(config.result_root_dir, "dr_{}.txt".format(data_id_0)),
        Y_0, labels)

    elapsed = 0
    # for each frame
    for i in range(1, len(data_ids)):
        data_id_1 = data_ids[i]
        data_dim_1 = data_dims[data_id_1]

        print("joint t-sne optimize {}".format(data_id_1))

        data_path = os.path.join(config.input_dir,
                                 "f_{}.txt".format(data_id_1))
        X1 = np.loadtxt(data_path, usecols=range(0, data_dim_1))
        labels = np.loadtxt(data_path, usecols=(data_dim_1, ), dtype=str)

        if config.norm:
            X1 = dr_utils.norm_data(X1)

        edges0 = dr_utils.read_graph(
            os.path.join(config.graph_dir, "g_{}.txt".format(data_id_0)),
            config.k_closest_count)
        edges1 = dr_utils.read_graph(
            os.path.join(config.graph_dir, "g_{}.txt".format(data_id_1)),
            config.k_closest_count)
        ''' discount point similarity by neighborhood ratio '''
        match_points = dr_utils.read_match_points(
            os.path.join(config.sim_dir,
                         "g_{}_g_{}.txt".format(data_id_0, data_id_1)))

        ratios = dr_utils.common_neighbor_ratio(
            node_num=X1.shape[0],
            edges0=edges0,
            edges1=edges1,
            Neighbor_size=config.k_closest_count)
        for vi in match_points:
            match_points[vi] *= ratios[vi]

        # compute edge similarity based on point similarity
        match_edges = dr_utils.find_match_edges(
            edges0=edges0, edges1=edges1, match_points=match_points)

        # normalize edges for better distinction
        edge_sim_min = np.min(list(match_edges.values()))
        edge_sim_max = np.max(list(match_edges.values()))
        if edge_sim_min != edge_sim_max:
            for e in match_edges:
                match_edges[e] = (match_edges[e] -
                                  edge_sim_min) / (edge_sim_max - edge_sim_min)


        start = time.perf_counter()
        '''note that we use the same initialization everytime'''
        Y1 = joint_tsne(Y_0=Y_0,
                        Y_1_I=Y_I,
                        X_1=X1,
                        match_edges=match_edges,
                        no_dims=2,
                        initial_dims_1=data_dim_1,
                        perplexity=perplexity,
                        gamma_=gamma,
                        verbose=0)  # perplexity

        elapsed += (time.perf_counter() - start)
        # save joint-tsne results
        dr_utils.saveVectors(
            os.path.join(config.result_root_dir, "dr_{}.txt".format(data_id_1)),
            Y1, labels)

        # the output is used to anchor next frame
        Y_0 = Y1
        data_id_0 = data_id_1

    print("joint t-sne optimization used:", elapsed)
    return elapsed


if __name__ == "__main__":
    argv = sys.argv
    assert (len(argv) == 2)
    config_path = argv[1]

    config_json = json.loads(open(config_path).read())
    config = SNEConfig(config_json)