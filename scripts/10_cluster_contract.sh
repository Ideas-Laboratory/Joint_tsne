# !/bin/bash
config_path="/Users/joe/Codes/PythonProjects/Joint_tsne/config/10_cluster_contract.json"

echo "configuration file: ${config_path}"

# 1. build knn graph
python3 ../codes/graphBuild/run.py $config_path

# 2. compute graph similarity
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
./graphSim.app/Contents/MacOS/graphSim $config_path
cd ../

# 3. run thesne
python3 ../codes/thesne/run.py $config_path
