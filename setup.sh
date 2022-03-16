#!/bin/bash -x

### get script directory
BIN=$(dirname "${BASH_SOURCE[0]}")
pushd $BIN > /dev/null
BASEDIR=`pwd`
popd > /dev/null

mkdir $BASEDIR/AF_data/
pushd $BASEDIR/AF_data/ > /dev/null

### download AF2 parameters
mkdir params
pushd params > /dev/null
wget https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar
tar -xf alphafold_params_2021-07-14.tar
rm alphafold_params_2021-07-14.tar
popd > /dev/null

### create singularity environment to run AF2
singularity pull AF_environment.sif docker://catgumag/alphafold:latest


# To get the latest from DeepMind you need to do this (as sudo)
#git clone git@github.com:/deepmind/alphafold.git
#cd alphafold
#sudo docker build -f docker/Dockerfile -t alphafold .
#sudo singularity build alphafold-multimer.sif docker-daemon://alphafold:latest
#popd > /dev/null


popd > /dev/null

### download  Uniclust30
wget http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/uniclust30_2018_08_hhsuite.tar.gz
tar -zxvf uniclust30_2018_08_hhsuite.tar.gz

pushd $BASEDIR/bin/ > /dev/null

### get DockQ program
git clone https://github.com/bjornwallner/DockQ.git
pushd DockQ/ > /dev/null
make
popd > /dev/null

### get CD-Hit program
git clone https://github.com/weizhongli/cdhit.git
pushd cdhit/ > /dev/null
make
popd > /dev/null

popd > /dev/null

