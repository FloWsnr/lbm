from pathlib import Path
import numpy as np
from skimage import measure
import scipy.ndimage.morphology


def extract_nw_fluid(structure: np.ndarray):
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


def remove_nw_fluid(structure: np.ndarray):
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


def get_surface(structure: np.ndarray):
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


def get_solid(structure: np.ndarray):
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

def add_inlet_outlet_layers(structure: np.ndarray, num_layers: int):
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

def create_geom_for_palabos(inputs):
    sim_dir = inputs["input output"]["simulation directory"]
    input_dir = inputs["input output"]["input folder"]
    geom_file_name = inputs["geometry"]["file name"]
    geom_file = sim_dir + "/" + input_dir + geom_file_name
    nx = inputs["domain"]["domain size"]["nx"]
    ny = inputs["domain"]["domain size"]["ny"]
    nz = inputs["domain"]["domain size"]["nz"]
    geom_name = inputs["domain"]["geom name"]

    # read-in file
    structure = np.load(geom_file)
    # select a subset for simulation
    structure = structure[0:nz, 0:ny, 0:nx]
    # Tanspose x and z
    structure = structure.transpose([2, 1, 0])

    # Create a mask for the non-wetting fluid
    nw_fluid_mask = extract_nw_fluid(structure)
    # Remove the non-wetting fluid from the structure
    structure = remove_nw_fluid(structure)

    # Create a surface geometry
    surface = get_surface(structure)
    # Create a solid geometry
    solid = get_solid(structure)

    # Turn all non-surface voxels into solid phase 2 (inner solid)
    structure[np.logical_and(solid, ~surface)] = 2

    # Add inlet and outlet layers
    num_layers = inputs["domain"]["inlet and outlet layers"]
    structure = add_inlet_outlet_layers(structure, num_layers)

    # Maybe structure must be converted to 2608 etc

    path = Path(sim_dir / input_dir)
    path.mkdir(parents=True, exist_ok=True)

    structure.flatten().tofile(path / f"{geom_name}.dat")  # Save geometry
    np.save(path / f"{geom_name}.npy", structure)
