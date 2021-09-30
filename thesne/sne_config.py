import os


# this class is irevelant from data generation
class SNEConfig:
    def __init__(self, config_json):
        # initlialize configurations
        # self.algo = AlgoConfig(config_json["algo"])
        self.perplexity = config_json["algo"]['perplexity']
        self.k_closest_count = config_json["algo"]['k_closest_count']
        #min(3*perplexity, pts_size)   # K近邻的个数+1（虽然是K=4，但由于包含自身，实际为K-1邻近）
        self.bfs_level = config_json["algo"]['bfs_level']
        self.gfd_calc = config_json["algo"]['gfd_calc']
        self.beta = config_json["algo"]['beta']
        self.gamma = config_json["algo"]['gamma']
        self.anchor_rate = config_json["algo"]['anchor_rate']
        self.dynamic_lambda = config_json["algo"]["dynamic_lambda"]

        self.data_name = config_json["thesne"]["data_name"]
        dirname = os.path.dirname(__file__)
        self.input_dir = os.path.join(dirname,
                                      "../../data/{}".format(self.data_name))
        self.graph_dir = os.path.join(
            dirname, "../../knn graph/{}".format(self.data_name))
        self.sim_dir = os.path.join(dirname,
                                    "../../graphSim/{}".format(self.data_name))

        method1 = "dynamic tsne"
        method2 = "equal initialization"
        method3 = "random initialization"
        method4 = "joint tsne"

        # output results directories
        self.output_root_dir = os.path.join(dirname, "../../results",
                                            self.data_name)
        self.dynamic_dir = os.path.join(self.output_root_dir, method1)
        self.equalinit_dir = os.path.join(self.output_root_dir, method2)
        self.randominit_dir = os.path.join(self.output_root_dir, method3)
        self.joint_tsne_dir = os.path.join(self.output_root_dir, method4)

        self.fig_root_dir = os.path.join(dirname, "../../figures",
                                         self.data_name)
        self.dynamic_fig_dir = os.path.join(self.fig_root_dir, method1)
        self.equalinit_fig_dir = os.path.join(self.fig_root_dir, method2)
        self.randominit_fig_dir = os.path.join(self.fig_root_dir, method3)
        self.joint_tsne_fig_dir = os.path.join(self.fig_root_dir, method4)

        # self.draw = DrawConfig(config_json["draw"])
        self.pts_size = config_json['thesne']['pts_size']
        self.data_ids = config_json['thesne']['data_ids']
        self.data_dims = config_json['thesne']['data_dims']
        self.data_titles = config_json['thesne']['data_titles']
        self.norm = config_json['thesne']['norm']
