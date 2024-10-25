from pathlib import Path
from typing import Optional, Union

import numpy as np
import scipy.ndimage.morphology
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
    shape = tuple(x - 1 for x in shape)

    loaded_data = {}
    for key in data.cell_data.keys():
        structure = data.cell_data[key]
        structure = structure.reshape(shape, order="F")
        loaded_data[key] = structure

    return loaded_data


class VTIWriter:
    def __init__(self, directory: Union[Path, str]) -> None:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        self.directory = directory

        self.vti_file = None
        self.grid = None

    def write(self, data: dict, file_name: str) -> None:
        vti_file = self.directory / file_name
        self.vti_file = vti_file
        if self.vti_file.exists():
            raise FileExistsError(f"File {self.vti_file} already exists")

        self.grid = pv.ImageData()
        for key, value in data.items():
            self._add_data(value, name=key)
        self.grid.save(self.vti_file)

    def _add_data(self, data: np.ndarray, name: str = "value") -> None:
        if data.ndim == 2:
            data = np.expand_dims(data, axis=2)
        shape = np.array(data.shape) + 1
        self.grid.dimensions = shape
        self.grid.origin = (0, 0, 0)
        self.grid.spacing = (1, 1, 1)
        self.grid.cell_data[name] = data.flatten(order="F")


class PalabosGeometry:
    def __init__(self, inputs: dict):
        sim_dir = inputs["input_output"]["simulation_directory"]
        sim_dir = Path(sim_dir)
        structure_file = sim_dir / inputs["input_output"]["file_name"]

        if structure_file.suffix == ".vti":
            data = read_vti_file(structure_file)
            org_structure = data["structure"]
        elif structure_file.suffix == ".npy":
            org_structure = np.load(structure_file)
        else:
            raise ValueError(f"Unknown file type: {structure_file.suffix}")

        # crop the structure
        self.structure = self._crop_structure(org_structure, inputs["geometry"])
        # swap x and z
        self.structure = np.swapaxes(self.structure, 0, 2)

        self.inout_layers = inputs["domain"]["inlet_outlet_layers"]
        self.materials = inputs["materials"]

        self.input_dir = inputs["input_output"]["input_folder"]
        self.output_dir = inputs["input_output"]["output_folder"]

    def _crop_structure(self, structure: np.ndarray, geom: dict) -> np.ndarray:
        """Crop the structure"""
        crop = geom["crop"]
        if crop:
            return structure[
                geom["x1"] : geom["x2"],
                geom["y1"] : geom["y2"],
                geom["z1"] : geom["z2"],
            ]
        else:
            return structure

    def _extract_nw_fluid(self, structure: np.ndarray) -> np.ndarray:
        """
        Extracts the non-wetting fluid from a structure

        Parameters
        ----------
        structure : np.ndarray
            The structure array

        Returns
        -------
        nw_fluid : np.ndarray
            The non-wetting fluid array

        """
        nw_fluid = np.zeros_like(structure)
        nw_fluid[structure == 3] = 3
        return nw_fluid

    def _remove_nw_fluid(self, structure: np.ndarray) -> np.ndarray:
        """
        Removes the non-wetting fluid from a structure

        Parameters
        ----------
        structure : np.ndarray
            The structure array

        Returns
        -------
        structure : np.ndarray
            The structure array with the non-wetting fluid removed

        """
        structure[structure == 3] = 0
        return structure

    def _get_surface(self, structure: np.ndarray) -> np.ndarray:
        """
        Gets the surface from a structure

        Parameters
        ----------
        structure : np.ndarray
            The structure array

        Returns
        -------
        surface : np.ndarray
            The surface array

        """

        distance_map: np.ndarray = scipy.ndimage.distance_transform_edt(
            structure, return_indices=False
        )
        surface = np.logical_and(distance_map > 0, distance_map < 2)
        return surface

    def _get_solid(self, structure: np.ndarray) -> np.ndarray:
        """
        Get a solid mask of the structure, i.e. all voxels that are solid are True

        Parameters
        ----------
        structure : np.ndarray
            The structure array

        Returns
        -------
        solid : np.ndarray
            The solid array

        """
        solid = np.zeros_like(structure, dtype=bool)
        solid[structure == 1] = True
        solid[structure == 2] = True
        solid[structure == 4] = True
        solid[structure == 6] = True
        solid[structure == 7] = True

        return solid

    def _add_boundary_bounces(self, structure: np.ndarray) -> np.ndarray:
        """Adds boundary bounces to a structure"""
        solid = self._get_solid(structure)

        structure[0, :, :] = 1
        structure[:, 0, :] = 1
        structure[:, :, 0] = 1

        structure[-1, :, :] = 1
        structure[:, -1, :] = 1
        structure[:, :, -1] = 1

        structure[~solid] = 0

        return structure

    def _add_inlet_outlet_layers(
        self, structure: np.ndarray, num_layers: int
    ) -> np.ndarray:
        """
        Adds inlet and outlet layers to a structure

        Parameters
        ----------
        structure : np.ndarray
            The structure array
        num_layers : int
            The number of layers to add

        Returns
        -------
        structure : np.ndarray
            The structure array with inlet and outlet layers added

        """
        structure = np.pad(
            structure,
            ((num_layers, num_layers), (0, 0), (0, 0)),
            mode="constant",
            constant_values=(0, 0),
        )
        return structure

    def convert_material_ids(self):
        """Convert the material ids (255, 200 etc) to the correct values for Palabos

        Parameters
        ----------
        conversion : dict, optional
            The conversion dictionary, by default None
        """

        for material in self.materials:
            old_id = material[0]
            new_id = material[1]
            self.structure[self.structure == old_id] = new_id

    def create_geom_for_palabos(self):
        # Create a mask for the non-wetting fluid
        nw_fluid_mask = self._extract_nw_fluid(self.structure)
        # Remove the non-wetting fluid from the structure
        self.structure = self._remove_nw_fluid(self.structure)

        # Create a surface geometry
        surface = self._get_surface(self.structure)
        # Create a solid geometry
        solid = self._get_solid(self.structure)

        # Turn all non-surface voxels into solid phase 2 (inner solid)
        self.structure[np.logical_and(solid, ~surface)] = 2

        # add boundary bounces backs
        self.structure = self._add_boundary_bounces(self.structure)

        # Add inlet and outlet layers
        num_layers = self.inout_layers
        self.structure = self._add_inlet_outlet_layers(self.structure, num_layers)

        self.input_dir.mkdir(parents=True, exist_ok=True)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        flat_structure = self.structure.flatten(order="C")
        np.savetxt(self.input_dir / "structure.dat", flat_structure, fmt="%d")

        np.save(self.input_dir / "structure.npy", self.structure)

        vti_writer = VTIWriter(self.input_dir)
        data = {"structure": self.structure}
        vti_writer.write(data, "structure.vti")
