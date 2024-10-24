import subprocess
from pathlib import Path
import yaml

import lbm_wetting.utils.structure_prep as prep
from lbm_wetting.utils.palabos_file import PalabosInputFile
import lbm_wetting.utils.pydantic_schemas as schemas


def parse_input_file(file: Path) -> dict:
    with open(file) as f:
        inputs = yaml.full_load(f)
    inputs = schemas.Config(**inputs).model_dump()
    inputs = _update_config(inputs)
    return inputs


def _update_config(config: dict) -> dict:
    sim_dir = config["input_output"]["simulation_directory"]
    sim_name = config["input_output"]["simulation_name"]

    config["input_output"]["input_folder"] = Path(sim_dir) / sim_name / "input"
    config["input_output"]["output_folder"] = Path(sim_dir) / sim_name / "output"
    return config


def run_2_phase_sim(config: dict):
    # Steps
    # 1) create geom for palabos
    # 2) create palabos input file
    # 3) run 2-phase sim

    # 1) Create Palabos geometry
    print("Creating efficient geometry for Palabos...")
    palabos_geom = prep.PalabosGeometry(config)
    palabos_geom.convert_material_ids()
    palabos_geom.create_geom_for_palabos()

    # 2) Create simulation input file
    print("Creating input file...")
    palabos_file = PalabosInputFile(config)
    palabos_file.create_input_file()

    # 3) Run 2-phase simulation
    print("Running 2-phase simulation...")
    num_procs = config["simulation"]["num_procs"]

    sim_directory = config["input_output"]["simulation_directory"]
    input_dir = config["input_output"]["input_folder"]
    input_folder = Path(sim_directory) / input_dir

    # Save config to input folder
    with open(input_folder / "config.yml", "w") as f:
        yaml.dump(config, f)

    input_file = input_folder / "2_phase_sim_input.xml"
    palabos_bin = Path(
        "/home/fw641779/Coding/lattice-boltzmann-wetting/mplbm-ut-mirror/src/2-phase_LBM/ShanChen"
    )

    subprocess.run(
        ["mpirun", "-np", f"{num_procs}", f"{palabos_bin}", f"{input_file}"], check=True
    )

    return


if __name__ == "__main__":
    input_config = Path(
        "/home/fw641779/Coding/lattice-boltzmann-wetting/lbm_wetting/twophase.yml"
    )
    config = parse_input_file(input_config)
    run_2_phase_sim(config)  # Run 2 phase sim
