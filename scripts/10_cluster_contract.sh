# !/bin/bash
config_path="config/10_cluster_contract.json"
echo "Configuration file: ${config_path}"


workdir=$(pwd)
echo "Working directory: ${workdir}"


# 1. build knn graph
python3 codes/graphBuild/run.py $config_path

# 2. compute graph similarity
buildDir="codes/graphSim/build"
if [ ! -d $buildDir ]; then
    mkdir $buildDir
    echo "create directory ${buildDir}"
else
    echo "directory ${buildDir} already exists."
fi
cd $buildDir
qmake ../
make

cd $workdir

# this only works in macOS, you can change it to binary name in your operating system
bin=$buildDir/graphSim.app/Contents/MacOS/graphSim
$bin $config_path

# 3. run thesne
python3 codes/thesne/run.py $config_path
