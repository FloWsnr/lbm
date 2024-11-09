"""
Prepares the 2-phase simulation.
This means converting the structure to a format that can be used by Palabos and creating the input file.

Creator: Florian Wiesner
Date: 2024-11-08
"""

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
    """
    Update the config with the simulation directory and name.
    If the simulation directory already exists, increment the sim_name.
    """
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
    """
    Prepare the 2-phase simulation.
    """
    print("  Creating efficient geometry for Palabos...")
    palabos_geom = prep.PalabosGeometry(config)
    palabos_geom.convert_material_ids()
    palabos_geom.create_geom_for_palabos()

    config["structure"] = palabos_geom.structure

    # 2) Create simulation input file
    print("  Creating input file...")
    palabos_file = PalabosInputFile(config)
    palabos_file.create_input_file()

    # Save config in input folder
    config_path = config["input_output"]["input_folder"] / "config.yml"
    _dump_config(config, config_path)


def _dump_config(config: dict, config_path: Path):
    """
    Dump the config to a file.
    """
    # Remove the structure from the config
    config.pop("structure")
    # convert paths to strings
    for key, value in config["input_output"].items():
        if isinstance(value, Path):
            config["input_output"][key] = str(value)

    with open(config_path, "w") as f:
        yaml.dump(config, f)


def _check_prep_already_exists(sim_dir: Path, sim_name: str) -> bool:
    """
    Check if the preparation already exists.
    """
    return (sim_dir / sim_name / "input" / "2_phase_sim_input.xml").exists()


# Defaults
config_path = Path(
    "/home/fw641779/Coding/lattice-boltzmann-wetting/lbm_wetting/twophase.yml"
)
sim_dir = Path("/hpcwork/fw641779/lbm/Toray-120C/55cov/structure0")
sim_name = "test_run_0"

parser = argparse.ArgumentParser()
parser.add_argument(
    "config",
    type=Path,
    nargs="?",
    default=config_path,
    help="Path to the configuration file",
)
parser.add_argument(
    "sim_dir",
    type=Path,
    nargs="?",
    default=sim_dir,
    help="Simulation directory",
)
parser.add_argument(
    "sim_name",
    type=str,
    nargs="?",
    default="test_run_0",
    help="Simulation name (default: test_run_0)",
)
args = parser.parse_args()

config_path = Path(args.config)
sim_dir = Path(args.sim_dir)
sim_name = args.sim_name

if _check_prep_already_exists(sim_dir, sim_name):
    print("Preparation already exists. Palabos simulation can start now.")
    exit(0)

print("Starting preparation of 2-phase simulation...")
config = parse_input_file(config_path)
config = update_config(config, sim_dir, sim_name)
prep_2_phase_sim(config)
print("Done preparing 2-phase simulation!")
print("----------------------------------")
