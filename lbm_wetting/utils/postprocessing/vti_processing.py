"""
Processes the vti files created by the simulation.
Reading and writing the vti files.
Writing also includes a pvd file for paraview to bundle the vti files.

Creator: Florian Wiesner
Date: 2024-11-08
"""

import numpy as np
from pathlib import Path
from typing import Union
import xml.etree.cElementTree as ET

import pyvista as pv
import yaml

from lbm_wetting.utils.postprocessing.video_writer import VideoWriter


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
        elif data.ndim == 3:
            shape = np.array(data.shape) + 1
        elif data.ndim == 4:
            shape = np.array(data.shape[:3]) + 1
        else:
            raise ValueError(f"Data has {data.ndim} dimensions, expected 2, 3 or 4")
        self.grid.dimensions = shape
        self.grid.origin = (0, 0, 0)
        self.grid.spacing = (1, 1, 1)

        if data.ndim == 2 or data.ndim == 3:
            self.grid.cell_data[name] = data.flatten(order="F")
        elif data.ndim == 4:
            self.grid.cell_data[name + "_x"] = data[:, :, :, 0].flatten(order="F")
            self.grid.cell_data[name + "_y"] = data[:, :, :, 1].flatten(order="F")
            self.grid.cell_data[name + "_z"] = data[:, :, :, 2].flatten(order="F")

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


class VTIProcessor:
    def __init__(self, config: dict) -> None:
        self.config = config
        self.sim_dir = Path(config["input_output"]["simulation_directory"])

        output_folder = config["input_output"]["output_folder"]
        self.output_folder = Path(output_folder)
        self.processed_folder = self.output_folder.parent / "processed"
        # create processed folder
        self.processed_folder.mkdir(parents=True, exist_ok=True)

    def _load_structure(self) -> np.ndarray:
        structure_file = self.sim_dir / "structure.vti"
        if not structure_file.exists():
            file_with_suffix = structure_file.with_suffix(".npy")
            structure = np.load(file_with_suffix)
        else:
            structure_data = read_vti_file(structure_file)
            structure = structure_data["structure"]
        return structure

    def _find_states(self) -> list[Path]:
        states = [f for f in self.output_folder.glob("rho_f1_*.vti")]
        # sort states first by number, i.e rho_f1_001_00002000.vti, 001 is the number
        # and then by time, i.e 00002000
        states.sort(
            key=lambda x: (int(x.stem.split("_")[2]), int(x.stem.split("_")[3]))
        )
        return states

    def _find_steady_states(self, states: list[Path]) -> list[int]:
        # steady states are the last state of each first component
        # i.e the last 001, last 002, etc. Find them by checking if the next file has a different first component
        steady_states = []
        for i, f in enumerate(states[:-1]):
            if states[i + 1].stem.split("_")[2] != f.stem.split("_")[2]:
                steady_states.append(i)
        steady_states.append(len(states) - 1)
        return steady_states

    def _find_velocity_files(self) -> list[Path]:
        velocity_files = [f for f in self.output_folder.glob("vel_f1_*.dat")]
        velocity_files.sort(key=lambda x: (int(x.stem.split("_")[2])))
        return velocity_files

    def process_vti_files(self) -> None:
        self.structure = self._load_structure()

        self.writer = VTIWriter(
            self.processed_folder, pvd_file=self.processed_folder / "sim.pvd"
        )

        self.steady_state_writer = VTIWriter(
            self.processed_folder, pvd_file=self.processed_folder / "steady_states.pvd"
        )

        self.states = self._find_states()
        self.steady_states = self._find_steady_states(self.states)

        for i, density_file in enumerate(self.states):
            comp1 = density_file.stem.split("_")[2]
            comp2 = density_file.stem.split("_")[3]
            file_name = "sim_" + comp1 + "_" + comp2 + ".vti"
            time_step = float(comp1 + "." + comp2)
            # check if file is already processed
            if (self.processed_folder / file_name).exists():
                print(f"File already processed: {density_file}")
                self.writer.write_pvd(file_name, time_step)
                continue
            print(f"Processing file: {density_file}")

            ref_structure = np.copy(self.structure)
            # Add inlet/output layers to the structure to match density shape
            num_layers = config["domain"]["inlet_outlet_layers"]
            ref_structure = np.pad(
                ref_structure,
                ((0, 0), (0, 0), (num_layers, num_layers)),
                mode="constant",
                constant_values=0,
            )

            density_data = read_vti_file(density_file)
            density = density_data["Density"]
            swap_xz = config["geometry"]["swap_xz"]
            if swap_xz:
                density = np.swapaxes(density, 0, 2)
            wet_structure = _convert_density(density, ref_structure, num_layers)

            data = {"wet_structure": wet_structure}
            self.writer.write_vti(data, file_name, time_step)

            if i in self.steady_states:
                self.steady_state_writer.write_pvd(file_name, time_step)

            # Delete the original density file
            density_file.unlink()

        self.writer.close()
        self.steady_state_writer.close()

    def process_velocity_files(self) -> None:
        self.velocity_writer = VTIWriter(
            self.processed_folder, pvd_file=self.processed_folder / "velocity.pvd"
        )

        self.structure = self._load_structure()
        ref_structure = np.copy(self.structure)
        # Add inlet/output layers to the structure to match density shape
        num_layers = config["domain"]["inlet_outlet_layers"]
        ref_structure = np.pad(
            ref_structure,
            ((0, 0), (0, 0), (num_layers, num_layers)),
            mode="constant",
            constant_values=0,
        )
        velocity_files = self._find_velocity_files()
        for file in velocity_files:
            print(f"Processing velocity file: {file}")
            data = np.loadtxt(file)
            # data = data.astype(np.float32)
            if self.config["geometry"]["swap_xz"]:
                x, y, z = ref_structure.shape
                data = data.reshape(z, y, x, 3)
                data = np.swapaxes(data, 0, 2)
            else:
                x, y, z = ref_structure.shape
                data = data.reshape(x, y, z, 3)

            data = {"velocity": data}
            timestep = int(file.stem.split("_")[2])
            self.velocity_writer.write_vti(data, file.stem + ".vti", timestep)

        # save the data
        self.velocity_writer.close()

    def write_video(self) -> None:
        video_writer = VideoWriter(
            self.processed_folder / "video.gif", fps=5, upsample_rate=2
        )

        states = [f for f in self.processed_folder.glob("sim_*.vti")]
        # sort states first by number, i.e sim_001_00002000.vti, 001 is the number
        # and then by time, i.e 00002000
        states.sort(
            key=lambda x: (int(x.stem.split("_")[1]), int(x.stem.split("_")[2]))
        )
        for file in states:
            data = read_vti_file(file)
            wet_structure = data["wet_structure"]
            slice_index = int(wet_structure.shape[0] / 2)
            wet_structure_slice = wet_structure[slice_index, :, :]
            video_writer.write(wet_structure_slice)

        video_writer.close()


def _convert_density(
    density: np.ndarray, structure: np.ndarray, num_layers: int
) -> np.ndarray:
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
    return config


if __name__ == "__main__":
    import argparse

    default_sim_dir = Path("/hpcwork/fw641779/lbm/Test_Cones")
    default_sim_name = "test_run_1"

    parser = argparse.ArgumentParser(
        description="Process VTI files for LBM wetting simulation."
    )
    parser.add_argument(
        "sim_dir",
        type=Path,
        nargs="?",
        default=default_sim_dir,
        help="Simulation directory",
    )
    parser.add_argument(
        "sim_name",
        type=str,
        nargs="?",
        default=default_sim_name,
        help="Simulation name",
    )
    args = parser.parse_args()

    sim_dir = Path(args.sim_dir)
    sim_name = args.sim_name

    print("Processing VTI files for 2-phase simulation...")

    config_path = sim_dir / sim_name / "input" / "config.yml"
    config = parse_input_file(config_path)

    processor = VTIProcessor(config)
    processor.process_vti_files()
    processor.process_velocity_files()
    processor.write_video()
    print("Processing completed.")
