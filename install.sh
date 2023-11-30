module purge
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load GCCcore/.11.3.0
module load CMake/3.24.3

cd ./MPLBM-UT/src
unzip palabos-v2.2.1.zip

cd 2-phase_LBM/build
cmake ..
make -j 2

cd ../../1-phase_LBM/build
cmake ..
make -j 2

# cd ../../../
# conda create -n mplbm python=3.10
# conda activate mplbm
# python -m pip install --upgrade pip
# python -m pip install --upgrade setuptools
# pip install ./MPLBM-UT/python/