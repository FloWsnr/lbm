import numpy as np
from pathlib import Path
from typing import Union
import xml.etree.cElementTree as ET

import pyvista as pv
import yaml

import lbm_wetting.utils.pydantic_schemas as schemas


def read_vti_file(path: Union[Path, str]) -> dict:
    """Reads a vti/vtk file and converts it to a cp.ndarray
    Only works for the vti files created by pyvista.

    Parameters
    ----------
    path : Union[Path, str]
        The path to the vti/vtk file.

    Returns
    -------
    loaded_data : dict
        A dictionary containing the data from the vti/vtk file
    """
    file = Path(path)
    data = pv.read(file)
    shape: tuple = data.dimensions
    cell_shape = tuple(x - 1 for x in shape)

    loaded_data = {}
    for key in data.cell_data.keys():
        structure = data.cell_data[key]
        structure = structure.reshape(cell_shape, order="F")
        loaded_data[key] = structure

    point_shape = tuple(x for x in shape)
    for key in data.point_data.keys():
        structure = data.point_data[key]
        structure = structure.reshape(point_shape, order="F")
        loaded_data[key] = structure

    return loaded_data


class VTIWriter:
    def __init__(self, directory: Union[Path, str], pvd_file: Path = None) -> None:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        self.directory = directory

        if pvd_file is not None:
            self.pvd_file = pvd_file
            self.pvd_root = ET.Element("VTKFile", type="Collection", version="0.1")
            self.pvd_collection = ET.SubElement(self.pvd_root, "Collection")
        else:
            self.pvd_file = None
            self.pvd_root = None
            self.pvd_collection = None

        self.grid = None

    def write_vti(self, data: dict, file_name: str, timestep: int = 0) -> None:
        vti_file = self.directory / file_name
        # make clean grid
        self.grid = pv.ImageData()
        for key, value in data.items():
            self._add_data(value, name=key)

        # save the grid
        self.grid.save(vti_file)
        self._save_timestep(vti_file, timestep)

    def _add_data(self, data: np.ndarray, name: str = "value") -> None:
        if data.ndim == 2:
            data = np.expand_dims(data, axis=2)
        shape = np.array(data.shape) + 1
        self.grid.dimensions = shape
        self.grid.origin = (0, 0, 0)
        self.grid.spacing = (1, 1, 1)
        self.grid.cell_data[name] = data.flatten(order="F")

    def write_pvd(self, file_name: str, timestep: int = 0) -> None:
        """Add a file to the pvd collection without writing a vti file"""
        vti_file = self.directory / file_name
        self._save_timestep(vti_file, timestep)

    def _save_timestep(self, vti_file: Path, timestep: int) -> None:
        if self.pvd_root is not None:
            timestep = ET.SubElement(
                self.pvd_collection,
                "DataSet",
                timestep=str(timestep),
                file=str(vti_file.relative_to(self.directory)),
            )

    def close(self) -> None:
        if self.pvd_root is not None:
            self._save_pvd()

    def _save_pvd(self) -> None:
        tree = ET.ElementTree(self.pvd_root)
        tree.write(self.pvd_file)


def process_vti_files(config: dict) -> None:
    sim_dir = Path(config["input_output"]["simulation_directory"])

    output_folder = config["input_output"]["output_folder"]
    processed_folder = output_folder.parent / "processed"
    processed_folder.mkdir(parents=True, exist_ok=True)

    assert processed_folder.exists(), "Processed folder does not exist"

    # check if structure.vti exists
    structure_file = sim_dir / "structure.vti"
    if not structure_file.exists():
        # try to find file with suffix .npy
        file_with_suffix = structure_file.with_suffix(".npy")
        if not file_with_suffix.exists():
            raise FileNotFoundError(f"Structure file not found: {structure_file}")
        structure_file = file_with_suffix
        structure = np.load(structure_file)
    else:
        structure_data = read_vti_file(structure_file)
        structure = structure_data["structure"]

    writer = VTIWriter(processed_folder, pvd_file=processed_folder / "sim.pvd")
    steady_state_writer = VTIWriter(
        processed_folder, pvd_file=processed_folder / "steady_states.pvd"
    )

    states = [f for f in output_folder.glob("rho_f1_*.vti")]
    # sort states first by number, i.e rho_f1_001_00002000.vti, 001 is the number
    # and then by time, i.e 00002000
    states.sort(key=lambda x: (int(x.stem.split("_")[2]), int(x.stem.split("_")[3])))

    # steady states are the last state of each first component
    # i.e the last 001, last 002, etc. Find them by checking if the next file has a different first component
    steady_states = []
    for i, f in enumerate(states[:-1]):
        if states[i + 1].stem.split("_")[2] != f.stem.split("_")[2]:
            steady_states.append(i)
    steady_states.append(len(states) - 1)

    for i, file in enumerate(states):
        density_data = read_vti_file(file)
        density = density_data["Density"]

        # slice density to remove the input/output layers from the z axes
        swap_xz = config["geometry"]["swap_xz"]
        num_layers = config["domain"]["inlet_outlet_layers"]
        density = density[num_layers:-num_layers, :, :]
        if swap_xz:
            density = np.swapaxes(density, 0, 2)
        wet_structure = _convert_density(density, structure)

        data = {"wet_structure": wet_structure}

        comp1 = file.stem.split("_")[2]
        comp2 = file.stem.split("_")[3]
        file_name = "sim_" + comp1 + "_" + comp2 + ".vti"
        time_step = float(comp1 + "." + comp2)
        writer.write_vti(data, file_name, time_step)

        if i in steady_states:
            steady_state_writer.write_pvd(file_name, time_step)

    writer.close()
    steady_state_writer.close()


def _convert_density(density: np.ndarray, structure: np.ndarray) -> np.ndarray:
    water = density > 1
    solid = structure != 0

    wet_structure = np.zeros_like(structure)
    wet_structure[solid] = structure[solid]
    wet_structure[water] = 128
    wet_structure = wet_structure.astype(np.uint8)

    return wet_structure


def parse_input_file(file: Path) -> dict:
    with open(file) as f:
        config = yaml.full_load(f)
    config = schemas.Config(**config).model_dump()
    return config


def update_config(config: dict, sim_dir: Path, sim_name: str) -> dict:
    config["input_output"]["simulation_directory"] = sim_dir
    config["input_output"]["simulation_name"] = sim_name
    config["input_output"]["input_folder"] = Path(sim_dir) / sim_name / "input"
    config["input_output"]["output_folder"] = Path(sim_dir) / sim_name / "output"

    return config


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Process VTI files for LBM wetting simulation."
    )
    parser.add_argument("config", type=Path, help="Path to the configuration file")
    parser.add_argument("sim_dir", type=Path, help="Simulation directory")
    parser.add_argument("sim_name", type=str, help="Simulation name")
    args = parser.parse_args()

    config_path = Path(args.config)
    sim_dir = Path(args.sim_dir)
    sim_name = args.sim_name

    print("Processing VTI files for 2-phase simulation...")
    print("  config:", config_path)
    print("  sim_dir:", sim_dir)
    print("  sim_name:", sim_name)

    config = parse_input_file(config_path)
    config = update_config(config, sim_dir, sim_name)

    process_vti_files(config)
