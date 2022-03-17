# setup w/ conda env alphafold for non-docker version AF2.0/2.1/2.2

env_file=/mnt/data/envs/.jpasrc
if [[ -f "$env_file" ]]; then exit 1; fi

external_repo_pth=/repo

set -e
conda activate alphafold
# search for needle tool required by DockQ
if [[ "$(which needle)" == "" ]]; then
    conda install -y emboss -c bioconda
fi

# cd-hit
pushd $external_repo_pth
git clone https://github.com/weizhongli/cdhit.git
cd cdhit && make
echo export PATH=$PWD':$PATH' >>$env_file
popd

# DockQ
pushd $external_repo_pth
git clone https://github.com/bjornwallner/DockQ.git
cd DockQ && make
echo export PATH=$PWD':$PATH' >>$env_file
popd

pushd $external_repo_pth
git clone https://github.com/YaoYinYing/FoldDock
