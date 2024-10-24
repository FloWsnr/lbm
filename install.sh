module purge
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load CMake/3.24.3

cd ./mplbm-ut-mirror/src
unzip palabos-v2.2.1.zip

cd 2-phase_LBM/build
cmake ..
make -j 2

cd ../../1-phase_LBM/build
cmake ..
make -j 2