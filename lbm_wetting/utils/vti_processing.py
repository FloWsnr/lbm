import numpy as np
from pathlib import Path
from typing import Union
import xml.etree.cElementTree as ET

import pyvista as pv


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

        self.vti_file = None
        self.grid = None

    def write(self, data: dict, file_name: str, timestep: int = 0) -> None:
        self._open_file(file_name)
        for key, value in data.items():
            self._add_data(value, name=key)
        self._save_timestep(timestep)

    def write_pvd(self, file_name: str, timestep: int = 0) -> None:
        """Add a file to the pvd collection without writing a vti file"""
        self._open_file(file_name)
        self._save_timestep(timestep)

    def close(self) -> None:
        if self.pvd_root is not None:
            self._save_pvd()

    def _open_file(self, file_name: str) -> None:
        vti_file = self.directory / file_name
        self.vti_file = vti_file

        # if self.vti_file.exists():
        #     raise FileExistsError(f"File {self.vti_file} already exists")

        self.grid = pv.ImageData()

    def _add_data(self, data: np.ndarray, name: str = "value") -> None:
        if data.ndim == 2:
            data = np.expand_dims(data, axis=2)
        shape = np.array(data.shape) + 1
        self.grid.dimensions = shape
        self.grid.origin = (0, 0, 0)
        self.grid.spacing = (1, 1, 1)
        self.grid.cell_data[name] = data.flatten(order="F")

    def _save_timestep(self, timestep: int) -> None:
        self.grid.save(self.vti_file)
        if self.pvd_root is not None:
            timestep = ET.SubElement(
                self.pvd_collection,
                "DataSet",
                timestep=str(timestep),
                file=str(self.vti_file.relative_to(self.directory)),
            )

    def _save_pvd(self) -> None:
        tree = ET.ElementTree(self.pvd_root)
        tree.write(self.pvd_file)


def process_vti_files(directory: Union[Path, str]) -> None:
    directory = Path(directory)
    output_folder = directory / "output"
    processed_folder = directory / "processed"
    processed_folder.mkdir(parents=True, exist_ok=True)

    structure_data = read_vti_file(directory.parent / "structure.vti")
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

        wet_structure = _convert_density(density, structure)

        data = {"wet_structure": wet_structure}

        comp1 = file.stem.split("_")[2]
        comp2 = file.stem.split("_")[3]
        file_name = "sim_" + comp1 + "_" + comp2 + ".vti"
        time_step = float(comp1 + "." + comp2)
        writer.write(data, file_name, time_step)

        if i in steady_states:
            steady_state_writer.write_pvd(file_name, time_step)

    writer.close()
    steady_state_writer.close()


def _convert_density(density: np.ndarray, structure: np.ndarray) -> np.ndarray:
    # transpose density to match the structure
    density = np.swapaxes(density, 0, 2)
    # remove the input/output layers from the z axes
    density = density[:, :, 4:-4]

    water = density > 1
    solid = structure != 0

    wet_structure = np.zeros_like(structure)
    wet_structure[solid] = structure[solid]
    wet_structure[water] = 128
    wet_structure = wet_structure.astype(np.uint8)

    return wet_structure


if __name__ == "__main__":
    directory = Path("/hpcwork/fw641779/lbm/Toray-120C/55cov/structure0/test_run")
    process_vti_files(directory)
