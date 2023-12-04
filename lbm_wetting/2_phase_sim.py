import subprocess
import numpy as np
import matplotlib.pyplot as plt
import mplbm_utils as mplbm
import structure_prep as prep


def run_2_phase_sim(inputs):
    # Steps
    # 1) create geom for palabos
    # 2) create palabos input file
    # 3) run 2-phase sim

    sim_directory = inputs["input output"]["simulation directory"]

    # 1) Create Palabos geometry
    print("Creating efficient geometry for Palabos...")
    palabos_geom = prep.PalabosGeometry(inputs)
    palabos_geom.convert_material_ids()
    palabos_geom.create_geom_for_palabos()

    # 2) Create simulation input file
    print("Creating input file...")
    mplbm.create_palabos_input_file(inputs)

    # 3) Run 2-phase simulation
    print("Running 2-phase simulation...")
    num_procs = inputs["simulation"]["num procs"]
    input_dir = inputs["input output"]["input folder"]
    simulation_command = f"mpirun -np {num_procs} ../../src/2-phase_LBM/ShanChen {input_dir}2_phase_sim_input.xml"
    file = open(f"{sim_directory}/{input_dir}run_shanchen_sim.sh", "w")
    file.write(f"{simulation_command}")
    file.close()

    simulation_command_subproc = f"bash {sim_directory}/{input_dir}run_shanchen_sim.sh"
    subprocess.run(simulation_command_subproc.split(" "))

    return


def run_rel_perm_sim(inputs):
    # Rel Perm Steps
    # 1) create 1-phase palabos input file for relperms
    # 2) create geoms for rel perm
    # 3) run 1-phase sim to get rel perms

    sim_directory = inputs["input output"]["simulation directory"]

    # 1) Create geoms for rel perm
    print("Creating rel perm geometries...")
    user_geom_name = inputs["domain"]["geom name"]
    inputs = mplbm.create_geom_for_rel_perm(inputs)

    # 2) Create simulation input file
    print("Creating input file...")
    inputs["simulation type"] = "rel perm"
    mplbm.create_palabos_input_file(inputs)
    inputs["domain"]["geom name"] = user_geom_name

    # 3) Write rel perm bash file
    print("Running rel perm simulation...")
    num_procs = inputs["simulation"]["num procs"]
    input_dir = inputs["input output"]["input folder"]
    output_dir = inputs["input output"]["output folder"]
    simulation_command = f"mpirun -np {num_procs} ../../src/1-phase_LBM/permeability {input_dir}relperm_input.xml"
    file = open(f"{sim_directory}/{input_dir}run_relperm_sim.sh", "w")
    file.write(f"{simulation_command}")
    file.close()

    simulation_command_subproc = f"bash {sim_directory}/{input_dir}run_relperm_sim.sh"
    make_4relperm_folder = f"mkdir {sim_directory}/{output_dir}/4relperm"
    subprocess.run(make_4relperm_folder.split(" "))
    subprocess.run(simulation_command_subproc.split(" "))

    return


def process_and_plot_results(inputs):
    # Process data
    mplbm.create_pressure_data_file(inputs)
    mplbm.create_relperm_data_file(inputs)

    # Load and plot data
    sim_dir = inputs["input output"]["simulation directory"] + "/"
    output_dir = inputs["input output"]["output folder"]
    Sw = np.loadtxt(f"{sim_dir + output_dir}data_Sw.txt")
    Pc = np.loadtxt(f"{sim_dir + output_dir}data_Pc.txt")
    krw = np.loadtxt(f"{sim_dir + output_dir}data_krw.txt")
    krnw = np.loadtxt(f"{sim_dir + output_dir}data_krnw.txt")

    plt.figure(figsize=[10, 8])
    mplbm.plot_capillary_pressure_data(Sw, Pc)
    plt.savefig(sim_dir + "pc_curve.png", dpi=300)

    plt.figure(figsize=[10, 8])
    mplbm.plot_rel_perm_data(Sw, krw, krnw)
    # plt.yscale('log')
    plt.savefig(sim_dir + "relperm_curve.png", dpi=300)

    plt.figure(figsize=[6, 10])
    mplbm.plot_pc_and_rel_perm(Sw, Pc, krw, krnw)
    plt.savefig(sim_dir + "pc_and_relperm_curve.png", dpi=300)

    return


if __name__ == "__main__":
    input_file = "/work/fw641779/wetting/structures/test/input/input.yml"
    inputs = mplbm.parse_input_file(input_file)  # Parse inputs
    run_2_phase_sim(inputs)  # Run 2 phase sim
    # run_rel_perm_sim(inputs)  # Run rel perm
    # process_and_plot_results(inputs)  # Plot results

    # plt.show()
