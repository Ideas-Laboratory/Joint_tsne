import os


# this class is irevelant from data generation
class SNEConfig:
    def __init__(self, config_json):
        # initlialize configurations
        self.perplexity = config_json["algo"]['perplexity']
        self.k_closest_count = config_json["algo"]['k_closest_count']
        self.bfs_level = config_json["algo"]['bfs_level']
        self.gamma = config_json["algo"]['gamma']

        self.data_name = config_json["thesne"]["data_name"]
        dirname = os.path.dirname(__file__)
        self.input_dir = os.path.join(dirname,
                                      "../../data/{}".format(self.data_name))
        self.graph_dir = os.path.join(
            dirname, "../../knn graph/{}".format(self.data_name))
        self.sim_dir = os.path.join(dirname,
                                    "../../graphSim/{}".format(self.data_name))

        # output results directories
        self.result_root_dir = os.path.join(dirname, "../../results",
                                            self.data_name)

        self.fig_root_dir = os.path.join(dirname, "../../figures",
                                         self.data_name)
        self.pts_size = config_json['thesne']['pts_size']
        self.data_ids = config_json['thesne']['data_ids']
        self.data_dims = config_json['thesne']['data_dims']
        self.data_titles = config_json['thesne']['data_titles']
        self.norm = config_json['thesne']['norm']
