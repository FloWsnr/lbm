from pathlib import Path
import numpy as np
from skimage import measure
import scipy.ndimage.morphology

class PalabosGeometry:

    def __init__(self, inputs):
        self.sim_dir = inputs["input output"]["simulation directory"]
        self.input_dir = inputs["input output"]["input folder"]

        self.geom_file_name = inputs["geometry"]["file name"]
        self.geom_file = Path(self.sim_dir) / self.input_dir  / self.geom_file_name
        self.geom_name = inputs["domain"]["geom name"]
        self.inout_layers = inputs["domain"]["inlet and outlet layers"]

        # read-in file
        self.structure: np.ndarray = np.load(self.geom_file)
        # Tanspose x and z
        self.structure = self.structure.transpose([2, 1, 0])

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

    def _add_inlet_outlet_layers(self, structure: np.ndarray, num_layers: int) -> np.ndarray:
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
        structure = np.pad(structure, [(num_layers, num_layers), (0, 0), (0, 0)],
                        mode="constant", constant_values=(0, 0))
        return structure


    def convert_material_ids(self, conversion: dict = None):
        """Convert the material ids (255, 200 etc) to the correct values for Palabos"""

        if conversion is None:
            conversion = {
                255: 1,
                200: 4,
                150: 6,
                128: 3, # This is the non-wetting fluid
            }

        for key, val in conversion.items():
            self.structure[self.structure == key] = val

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

        # Add inlet and outlet layers
        num_layers = self.inout_layers
        self.structure = self._add_inlet_outlet_layers(self.structure, num_layers)

        # Maybe structure must be converted to 2608 etc

        path = Path(self.sim_dir) / self.input_dir
        path.mkdir(parents=True, exist_ok=True)

        self.structure.flatten().tofile(path / f"{self.geom_name}.dat")  # Save geometry
        np.save(path / f"{self.geom_name}.npy", self.structure)