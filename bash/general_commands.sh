################## common terminal ###############
# release linux cache memory
sudo sh -c 'echo 1 >  /proc/sys/vm/drop_caches'

# check if GPU is okay
python -c "import torch; print(torch.cuda.is_available())"

# count how much files are under current directory
find . -maxdepth 1 -type d -print0 | while IFS= read -r -d '' i ; do
    echo -n $i": " ; 
    (find "$i" -type f | wc -l) ; 
done

# replace all dot with dash of filenames under a directory
for f in *.*; do pre="${f%.*}"; suf="${f##*.}"; mv -i -- "$f" "${pre//./_}.${suf}"; done

# connect to local colab
sudo docker run --rm --runtime=nvidia --gpus all ubuntu nvidia-smi


################## pymol ###############
# pymol reindex
alter (fold_2025_05_19_22_11_atrkd_chk1_sq2_atp_mg_model_0 and chain A),resi=str(int(resi)+2295)
# pymol change chain id
alter (chain B), chain='A'
# pymol align
align mobile, reference
align mobile and chain A, reference and chain B
# remove res
remove sol
remove resn NA,CL
hide hrydrogen
H_add


################## WSL ###############
# mount a windows drive to wsl
sudo mkdir /mnt/driveX
sudo mount -t drvfs M: /mnt/driveX

# mount a windows share to wsl if already have access to share
# use single quotes otherwise you will need to escape the backslashes
sudo mkdir /mnt/mint
sudo mount -t drvfs '\\10.36.172.157\mint' /mnt/mint

# remove a mounted driver
sudo rm /mnt/driveX


################## miniforge ##################
# install mamba

# export and import env
mamba env export --no-builds | grep -v "^prefix: " > ngs.yml
mamba env create -f ngs.yml -p ~/miniforge3/envs/ngs

# mamba export all envs
for env in $(mamba env list | cut -d" " -f1); do 
   if [[ ${env:0:1} == "#" ]] ; then continue; fi;
   mamba env export -n $env > ./${env}.yml
done

# mamba remove all installed package in base env
mamba remove `conda list|awk 'NR>3 {print $1}'|tr '\n' ' '`
mamba remove `conda list|awk 'NR>3 {print $1}'|tr '\n' ' '`

    
################## run shinny APP ##################
# run shinyAPP from terminal
cd ~/Documents/GitHub/myScript/R_scripts/shiny_generate_well_positions
nohup R -e "shiny::runApp('shinyApp.R')" &


################## local PC ##################
Anydesk Shulab_Desktop: 1936728718 123456abc-


################## samba share ##################
# start samba server
sudo service smbd restart
sudo ufw allow samba
# list users
sudo pdbedit -L -v 
# add samba user
sudo smbpasswd -a [username]
# change pwd
sudo smbpasswd [username]
# remove samba user
sudo pdbedit -x -u [username]


################## cellprofiler-analyst ##################
# run cellprofiler-analyst
conda activate cellprofiler-analyst
python3 -m CellProfiler-Analyst

# modify sqlite path manually for different OS usage, run below code in SQL supportive software, like sqlite db browser
UPDATE Per_Image SET Image_PathName_ch1="Z:/Data/Project_condensate_reporter/2025-03-30_JY/Images";
UPDATE Per_Image SET Image_PathName_ch2="Z:/Data/Project_condensate_reporter/2025-03-30_JY/Images";
UPDATE Per_Image SET Image_PathName_ch3="Z:/Data/Project_condensate_reporter/2025-03-30_JY/Images";
UPDATE Per_Image SET Image_PathName_ch4="Z:/Data/Project_condensate_reporter/2025-03-30_JY/Images";
UPDATE Per_Image SET Image_ObjectsPathName_mask_cp_masks_cell="Z:/Data/Project_condensate_reporter/2025-03-30_JY/Images";






################## install amber24 ##################
# ubuntu 22.04 lts, failed....

# require a conda env
# if no amber MD needed, only install ambertools24 with conda for rapid installment

# need cuda version <=12.4 for amber24
 
# install amber24/ambertools24 from source with GPU acceleration
# download from website, put it at ~ 
# ambertools24.tar.bz2 must be located at same dir as amber24.tar.bz2
# tar files will generate ambertools24_src directory (for building)
cd ~
tar xvf Amber24.tar.bz2 
tar xvf AmberTools24.tar.bz2

sudo apt update
sudo apt upgrade
sudo apt install -y \
   tcsh make gcc g++ gfortran cmake flex bison \
   patch bc wget xorg-dev libz-dev libbz2-dev 

conda create -n amber24 -y
conda activate amber24
conda install -y \
   numpy scipy cython ipython notebook matplotlib \
   setuptools
conda install -y \
   nvidia/label/cuda-12.4.0::cuda-toolkit \
   nvidia/label/cuda-12.4.0::cuda-nvcc
# if need multiple gpu or cluster or no gpu, install mpi lib
# conda install -y openmpi mpi4py


# make file: ambertools24_src/build/run_cmake
# replace part: # Assume this is Linux
# this is for using ourself env to build amber (single GPU support, no mpi)
  cmake $AMBER_PREFIX/amber24_src \
    -DCMAKE_INSTALL_PREFIX=$AMBER_PREFIX/amber24 \
    -DCOMPILER=GNU  \
    -DMPI=FALSE \
    -DCUDA=TRUE \
    -DINSTALL_TESTS=TRUE \
    -DDOWNLOAD_MINICONDA=FALSE -DUSE_CONDA_LIBS=TRUE \
    -DCUDA_TOOLKIT_ROOT_DIR=$HOME/miniconda3/envs/amber24 \
    -DPYTHON_INCLUDE_DIR=$(python3 -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())")  \
    -DPYTHON_LIBRARY=$(python3 -c "import distutils.sysconfig as sysconfig; print(sysconfig.get_config_var('LIBDIR'))") \
    2>&1 | tee  cmake.log


# dir fir placing final files
mkdir $HOME/amber24

# cmake 
cd amber24_src/build

# if error, clean and then re-do cmake, this is very important
./clean_build -y & ./run_cmake

# check cmake log, if no error, proceed to compile
# -j 24 multi cores need sudo permission 
sudo make install -j 24


# using boost compiler, will need add intermediate PATH
export '$HOME/amber24_src/AmberTools/src/boost:$HOME/amber24_src/build/AmberTools/src/boost/stage/lib:$PATH'

# check
source $HOME/amber24/amber.sh
sander --version
pmemd.cuda --version


# optional, add amber to path
# add source $HOME/amber24/amber.sh to .bashrc

# Q&A
# if some library is missing, it may need to re=do ./run_cmake after installing missing packages
# if compile can't find some packages but it indeed installed correctly, just copy to that build file directory 

PATH=$PATH:/usr/include & export PATH





