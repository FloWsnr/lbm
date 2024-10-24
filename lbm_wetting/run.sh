#!/usr/bin/zsh
module purge
module load GCC/11.3.0
module load OpenMPI/4.1.4

palabos_dir="/home/fw641779/Coding/lattice-boltzmann-wetting/mplbm-ut-mirror/src/2-phase_LBM/ShanChen"

config_file="/home/fw641779/Coding/lattice-boltzmann-wetting/lbm_wetting/twophase.yml"
sim_dir="/hpcwork/fw641779/lbm/Toray-120C/55cov/structure0"
sim_name="test_run"

conda activate lbm
python3 /home/fw641779/Coding/lattice-boltzmann-wetting/lbm_wetting/2_phase_sim.py $config_file $sim_dir $sim_name
return_code=$?
if [ $return_code -nq 0 ]; then
    echo "Python script failed with code $return_code. Exiting..."
    return $return_code
else
    # Else continue
    input_file=$sim_dir/$sim_name/input/2_phase_sim_input.xml
    $MPIEXEC $FLAGS_MPI_BATCH $palabos_dir $input_file
fi
