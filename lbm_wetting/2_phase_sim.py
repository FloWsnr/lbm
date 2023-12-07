import subprocess
import mplbm_utils as mplbm
from pathlib import Path
import structure_prep as prep


def run_2_phase_sim(inputs):
    # Steps
    # 1) create geom for palabos
    # 2) create palabos input file
    # 3) run 2-phase sim

    sim_directory = inputs["input output"]["simulation directory"]
    input_dir = inputs["input output"]["input folder"]

    # 1) Create Palabos geometry
    print("Creating efficient geometry for Palabos...")
    palabos_geom = prep.PalabosGeometry(inputs)
    palabos_geom.convert_material_ids()
    palabos_geom.create_geom_for_palabos()

    # 2) Create simulation input file
    print("Creating input file...")
    palabos_file = prep.PalabosInputFile(inputs)
    palabos_file.create_input_file()

    # 3) Run 2-phase simulation
    print("Running 2-phase simulation...")
    num_procs = inputs["simulation"]["num procs"]

    input_folder = Path(sim_directory) / input_dir
    input_file = input_folder / "2_phase_sim_input.xml"

    palabos_bin = Path(
        "/home/fw641779/Coding/lattice-boltzmann-wetting/mplbm-ut-mirror/src/2-phase_LBM/ShanChen"
    )

    subprocess.run(
        ["mpirun", "-np", f"{num_procs}", f"{palabos_bin}", f"{input_file}"], check=True
    )

    return


if __name__ == "__main__":
    input_file = "/work/fw641779/wetting/structures/test/input/input.yml"
    inputs = mplbm.parse_input_file(input_file)  # Parse inputs
    run_2_phase_sim(inputs)  # Run 2 phase sim
    # run_rel_perm_sim(inputs)  # Run rel perm
    # process_and_plot_results(inputs)  # Plot results

    # plt.show()
