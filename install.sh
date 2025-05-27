module purge
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load CMake/3.24.3

cd ./palabos_sim/src
unzip palabos-v2.2.1.zip

mkdir -p 2-phase_LBM/build
cd 2-phase_LBM/build
cmake ..
make -j 2

cd ../../1-phase_LBM
mkdir -p build
cd build
cmake ..
make -j 2