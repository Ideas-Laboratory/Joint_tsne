import numpy as np
import pylab
from matplotlib.colors import ListedColormap


if __name__ == "__main__":
    data_name = "10_cluster_contract"
    frame_index = 9

    input_path = "../../results/{}/dr_{}.txt".format(
        data_name, frame_index)

    Y = np.loadtxt(input_path, usecols=range(0, 2), encoding='utf-8')
    labels = np.loadtxt(input_path, usecols=(2, ))


    colormap = [
        "#aecde1", "#3b77af", "#bbdd93", "#559d3f", "#ee9f9c", "#d1352b",
        "#c6b3d4", "#644195", "#ffffa6", "#a65e34"
    ]
    my_cmap = ListedColormap(colormap)

    pylab.scatter(Y[:, 0], Y[:, 1], 30, labels, cmap=my_cmap)
    pylab.show()