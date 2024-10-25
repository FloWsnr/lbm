from pathlib import Path
import yaml
import argparse

import lbm_wetting.utils.structure_prep as prep
from lbm_wetting.utils.palabos_file import PalabosInputFile
import lbm_wetting.utils.pydantic_schemas as schemas


def parse_input_file(file: Path) -> dict:
    with open(file) as f:
        config = yaml.full_load(f)
    config = schemas.Config(**config).model_dump()
    return config


def update_config(config: dict, sim_dir: Path, sim_name: str) -> dict:
    sim_run_dir = Path(sim_dir) / sim_name
    if sim_run_dir.exists():
        print("  Simulation directory already exists. Incrementing sim_name.")
        sim_name = f"{sim_name}_1"

    config["input_output"]["simulation_directory"] = sim_dir
    config["input_output"]["simulation_name"] = sim_name
    config["input_output"]["input_folder"] = Path(sim_dir) / sim_name / "input"
    config["input_output"]["output_folder"] = Path(sim_dir) / sim_name / "output"

    return config


def prep_2_phase_sim(config: dict):
    print("  Creating efficient geometry for Palabos...")
    palabos_geom = prep.PalabosGeometry(config)
    palabos_geom.convert_material_ids()
    palabos_geom.create_geom_for_palabos()

    config["structure"] = palabos_geom.structure

    # 2) Create simulation input file
    print("  Creating input file...")
    palabos_file = PalabosInputFile(config)
    palabos_file.create_input_file()


parser = argparse.ArgumentParser()
parser.add_argument("config", type=Path, help="Path to the configuration file")
parser.add_argument("sim_dir", type=Path, help="Simulation directory")
parser.add_argument("sim_name", type=str, help="Simulation name")
args = parser.parse_args()

config_path = Path(args.config)
sim_dir = Path(args.sim_dir)
sim_name = args.sim_name

print("Starting preparation of 2-phase simulation...")
print("  config:", config_path)
print("  sim_dir:", sim_dir)
print("  sim_name:", sim_name)

config = parse_input_file(config_path)
config = update_config(config, sim_dir, sim_name)
prep_2_phase_sim(config)
print("Done preparing 2-phase simulation!")
print("----------------------------------")
