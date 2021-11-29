# Joint t-sne
This is the implementation for paper [Joint t-SNE for Comparable Projections of Multiple High-Dimensional Datasets](http://www.yunhaiwang.net/Vis2021/joint-tsne). 

## abstract:
We present Joint t-Stochastic Neighbor Embedding (Joint t-SNE), a technique to generate comparable projections of multiple high-dimensional datasets. Although t-SNE has been widely employed to visualize high-dimensional datasets from various domains, it is limited to projecting a single dataset. When a series of high-dimensional datasets, such as datasets changing over time, is projected independently using t-SNE, misaligned layouts are obtained. Even items with identical features across datasets are projected to different locations, making the technique unsuitable for comparison tasks. To tackle this problem, we introduce edge similarity, which captures the similarities between two adjacent time frames based on the Graphlet Frequency Distribution (GFD). We then integrate a novel loss term into the t-SNE loss function, which we call vector constraints, to preserve the vectors between projected points across the projections, allowing these points to serve as visual landmarks for direct comparisons between projections. Using synthetic datasets whose ground-truth structures are known, we show that Joint t-SNE outperforms existing techniques, including Dynamic t-SNE, in terms of local coherence error, Kullback-Leibler divergence, and neighborhood preservation. We also showcase a real-world use case to visualize and compare the activation of different layers of a neural network.


### Environment
+ This is a hybrid programming based on C++ and Python, and supported by shell script.
+ It requires [Qt](https://www.qt.io/), [Python 3.6](https://www.python.org/), [numpy](https://numpy.org/) and [scikit-learn](https://scikit-learn.org/).

## How to use
1. Put the directory of your data sequence, e.g. "YOUR_DATA" in **Joint_tsne/data**. There are several requirements on the format and organization of your data: 
   + Each data frame is named as *f_i.txt*, where *i* is the time step/index of this data frame in the sequence.
   + The *j* th row of the data frame contains both the feature vector and label of the *j* th item, which is seperated by \tab. The label is at the last position.
   + All data frames must have the same number of rows, and the the same item is at the same row in different data frames to compute the node similarities one by one.  


2. Create a configuration file, e.g. "YOUR_DATA.json" in **Joint_tsne/config**, which is organized as a json structure.

<code>

```json
{
  "algo": {
    "k_closest_count": 3,
    "perplexity": 70,
    "bfs_level": 1,
    "gamma": 0.1
  },
  "thesne": {
    "data_name": "YOUR_DATA",
    "pts_size": 2000,
    "norm": false,
    "data_ids": [1, 3, 6, 9],
    "data_dims": [100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
    "data_titles": [
      "t=0",
      "t=1",
      "t=2",
      "t=3",
      "t=4",
      "t=5",
      "t=6",
      "t=7",
      "t=8",
      "t=9"
    ]
  }
}
```
</code>

In this file, *algo* represents the hyperparamters of our algorithm except for *bfs_level*, which always equals to 1. *thesne* contains the information of the input data. Please remember that *data_name* must be consistent with the directory name in the previous step.

3. Create a shell script, e.g. "YOUR_DATA.sh" in **Joint_tsne/scripts** as below:

<code>
    
```shell
# !/bin/bash
# 1. specify the configuration file with absolute file path
config_path="xxx/Joint_tsne/config/YOUR_DATA.json"

workdir=$(cd $(dirname $0); pwd)

# 2. build knn graph for each data frame
python3 ../codes/graphBuild/run.py $config_path

# 3. compute edge similarities between each two adjacent data frames
buildDir="../codes/graphSim/build"
if [ ! -d $buildDir ]; then
    mkdir $buildDir
    echo "create directory ${buildDir}"
else
    echo "directory ${buildDir} already exists."
fi
cd $buildDir
qmake ../
make

# bin is dependent on your operating system
bin=./graphSim.app/Contents/MacOS/graphSim
$bin $config_path

cd $workdir

# 4. run t-sne optimization
python3 ../codes/thesne/run.py $config_path
```
</code>

There are several places you should pay attention to. 
+ Again, *config_path* must be consitent with the name of configuration file in previous step
+ *bin* is dependent on your operating system. If you use linux, you should change it to 

        bin=./graphSim

4. change your directory to **Joint_tsne/scripts** and type 

<code>

    sh YOUR_DATA.sh

</code>

The final embeddings will be generated in **Joint_tsne/results/YOUR_DATA**.




### Example
You can find an example in **Joint_tsne/scripts/10_cluster_contract.sh**.