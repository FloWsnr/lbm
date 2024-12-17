from pathlib import Path

import numpy as np
import scipy.ndimage.morphology

from lbm_wetting.utils.postprocessing.vti_processing import VTIWriter, read_vti_file


def _load_structure(structure_file: Path) -> np.ndarray:
    if structure_file.suffix == ".vti":
        data = read_vti_file(structure_file)
        org_structure = data["structure"]
    elif structure_file.suffix == ".npy":
        org_structure = np.load(structure_file)
    else:
        raise FileNotFoundError(f"Structure file not found: {structure_file}")
    return org_structure


class PalabosGeometry:
    def __init__(self, inputs: dict):
        sim_dir = inputs["input_output"]["simulation_directory"]
        sim_dir = Path(sim_dir)

        # Search for the structure file
        file_name = inputs["input_output"]["file_name"]
        files = list(sim_dir.glob(f"{file_name}.*"))
        if len(files) == 0:
            raise FileNotFoundError(f"Structure file not found: {file_name}")
        elif len(files) > 1:
            raise FileNotFoundError(f"Multiple structure files found: {file_name}")
        structure_file = files[0]

        org_structure = _load_structure(structure_file)

        # crop the structure
        self.structure = self._crop_structure(org_structure, inputs["geometry"])

        if inputs["geometry"]["swap_xz"]:
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
        Is only needed if you want water present in the structure at the beginning of the simulation.

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
        Gets the surface from a structure.
        This means all voxels that are solid but directly next to the pores.

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
        """Adds boundary bounces to a structure.
        Places a layer of solid material (phase 1) around the structure.
        Is needed so that all solid voxels at the surface are solid phase 1 and not phase 2.
        Phase 2 is the inner solid phase.
        """
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
        self, structure: np.ndarray, num_layers: int, single_phase: bool = False
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
        if single_phase:
            structure = np.pad(
                structure,
                ((num_layers, num_layers), (0, 0), (0, 0)),
                mode="constant",
                constant_values=(0, 0),
            )
        else:
            structure = np.pad(
                structure,
                ((num_layers, num_layers), (0, 0), (0, 0)),
                mode="constant",
                constant_values=(3, 0),
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

    def create_geom_for_palabos(self, single_phase: bool = False):
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
        self.structure = self._add_inlet_outlet_layers(
            self.structure, num_layers, single_phase
        )

        self.input_dir.mkdir(parents=True, exist_ok=True)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        flat_structure = self.structure.flatten(order="C")
        np.savetxt(self.input_dir / "structure.dat", flat_structure, fmt="%d")

        vti_writer = VTIWriter(self.input_dir)
        data = {"structure": self.structure}
        vti_writer.write_vti(data, "structure.vti")
