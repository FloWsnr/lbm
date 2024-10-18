from pathlib import Path
from typing import Optional
import numpy as np
import scipy.ndimage.morphology


class PalabosGeometry:
    def __init__(self, inputs: dict):
        self.sim_dir = inputs["input_output"]["simulation_directory"]
        self.input_dir = inputs["input_output"]["input_folder"]

        self.geom_file_name = inputs["geometry"]["file_name"]
        self.geom_file = Path(self.sim_dir) / self.input_dir / self.geom_file_name
        self.geom_name = inputs["domain"]["geom_name"]
        self.inout_layers = inputs["domain"]["inlet_outlet_layers"]

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

    def convert_material_ids(self, conversion: Optional[dict] = None):
        """Convert the material ids (255, 200 etc) to the correct values for Palabos

        Parameters
        ----------
        conversion : dict, optional
            The conversion dictionary, by default None
        """

        if conversion is None:
            conversion = {
                255: 1,
                200: 4,
                150: 6,
                128: 3,  # This is the non-wetting fluid
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

        path: Path = Path(self.sim_dir) / self.input_dir
        path.mkdir(parents=True, exist_ok=True)

        self.structure.flatten().tofile(path / f"{self.geom_name}.dat")  # Save geometry
        np.save(path / f"{self.geom_name}.npy", self.structure)


class PalabosInputFile:
    def __init__(self, inputs):
        self.inputs = inputs

        sim_dir = inputs["input_output"]["simulation_directory"]
        input_dir = inputs["input_output"]["input_folder"]
        output_folder = inputs["input_output"]["output_folder"]

        self.output_dir = Path(sim_dir) / output_folder

        self.file = Path(sim_dir) / input_dir / "2_phase_sim_input.xml"

        restart_sim = inputs["simulation"]["restart_sim"]
        with open(self.file, "w") as file:
            file.write('<?xml version="1.0"?>\n\n')
            file.write(f"<load_savedstated> {restart_sim} </load_savedstated>\n\n")

    def _write_geometry_section(self):
        # Get geometry inputs
        periodic_x = self.inputs["domain"]["periodic_boundary"]["x"]
        periodic_y = self.inputs["domain"]["periodic_boundary"]["y"]
        periodic_z = self.inputs["domain"]["periodic_boundary"]["z"]

        nx = self.inputs["domain"]["domain_size"]["nx"]
        ny = self.inputs["domain"]["domain_size"]["ny"]
        nz = self.inputs["domain"]["domain_size"]["nz"]
        num_layers = self.inputs["domain"]["inlet_outlet_layers"]
        domain_size = [nx + (2 * num_layers), ny, nz]

        geo_file_name = self.inputs["domain"]["geom_name"]
        geo_file = self.file.parent / f"{geo_file_name}"

        with open(self.file, "a") as file:
            # Write geometry section
            file.write("<geometry>\n")
            # Geometry name
            file.write(f"\t<file_geom> {geo_file} </file_geom>\n")
            # Geometry size
            file.write(
                f"\t<size> <x> {domain_size[0]} </x> <y> {domain_size[1]} </y> <z> {domain_size[2]} </z> </size>\n"
            )
            # Periodicity
            file.write(f"\t<per>\n")
            file.write(
                f"\t\t<fluid1> <x> {periodic_x} </x> <y> {periodic_y} </y> <z> {periodic_z} </z> </fluid1>\n"
            )
            file.write(
                f"\t\t<fluid2> <x> {periodic_x} </x> <y> {periodic_y} </y> <z> {periodic_z} </z> </fluid2>\n"
            )
            file.write(f"\t</per>\n")
            file.write("</geometry>\n\n")

    def _write_fluid_positions(self):
        load_fluid_type = self.inputs["simulation"]["fluid_init"]
        if load_fluid_type == "geom":
            load_fluid_from_geom = True
        else:
            load_fluid_from_geom = False

        fluid1_x1 = self.inputs["simulation"]["fluid_1_init"]["x1"]
        fluid1_x2 = self.inputs["simulation"]["fluid_1_init"]["x2"]
        fluid1_y1 = self.inputs["simulation"]["fluid_1_init"]["y1"]
        fluid1_y2 = self.inputs["simulation"]["fluid_1_init"]["y2"]
        fluid1_z1 = self.inputs["simulation"]["fluid_1_init"]["z1"]
        fluid1_z2 = self.inputs["simulation"]["fluid_1_init"]["z2"]
        fluid2_x1 = self.inputs["simulation"]["fluid_2_init"]["x1"]
        fluid2_x2 = self.inputs["simulation"]["fluid_2_init"]["x2"]
        fluid2_y1 = self.inputs["simulation"]["fluid_2_init"]["y1"]
        fluid2_y2 = self.inputs["simulation"]["fluid_2_init"]["y2"]
        fluid2_z1 = self.inputs["simulation"]["fluid_2_init"]["z1"]
        fluid2_z2 = self.inputs["simulation"]["fluid_2_init"]["z2"]

        with open(self.file, "a") as file:
            # Write initial position of fluids
            file.write(f"<init>\n")
            file.write(
                f"\t<fluid_from_geom> {load_fluid_from_geom} </fluid_from_geom>\n"
            )
            file.write(f"\t<fluid1>\n")
            file.write(
                f"\t\t <x1> {fluid1_x1} </x1> <y1> {fluid1_y1} </y1> <z1> {fluid1_z1} </z1>\n"
            )
            file.write(
                f"\t\t <x2> {fluid1_x2} </x2> <y2> {fluid1_y2} </y2> <z2> {fluid1_z2} </z2>\n"
            )
            file.write(f"\t</fluid1>\n")
            file.write(f"\t<fluid2>\n")
            file.write(
                f"\t\t <x1> {fluid2_x1} </x1> <y1> {fluid2_y1} </y1> <z1> {fluid2_z1} </z1>\n"
            )
            file.write(
                f"\t\t <x2> {fluid2_x2} </x2> <y2> {fluid2_y2} </y2> <z2> {fluid2_z2} </z2>\n"
            )
            file.write(f"\t</fluid2>\n")
            file.write("</init>\n\n")

    def _write_fluid_data(self):
        rho_f1 = self.inputs["simulation"]["rho_f1"]
        rho_f2 = self.inputs["simulation"]["rho_f2"]

        force_f1 = self.inputs["simulation"]["force_f1"]
        force_f2 = self.inputs["simulation"]["force_f2"]

        Gc = self.inputs["simulation"]["fluid_data"]["Gc"]
        omega_f1 = self.inputs["simulation"]["fluid_data"]["omega_f1"]
        omega_f2 = self.inputs["simulation"]["fluid_data"]["omega_f2"]
        G_ads_f1_s1 = self.inputs["simulation"]["fluid_data"]["G_ads_f1_s1"]
        G_ads_f1_s2 = self.inputs["simulation"]["fluid_data"]["G_ads_f1_s2"]
        G_ads_f1_s3 = self.inputs["simulation"]["fluid_data"]["G_ads_f1_s3"]
        G_ads_f1_s4 = self.inputs["simulation"]["fluid_data"]["G_ads_f1_s4"]

        pressure_bc = self.inputs["simulation"]["pressure_bc"]
        if pressure_bc == True:
            minimum_radius = self.inputs["simulation"]["minimum_radius"]
            num_pc_steps = self.inputs["simulation"]["num_pressure_steps"]
        else:
            minimum_radius = 1
            num_pc_steps = 0

        with open(self.file, "a") as file:
            # Write fluid data
            file.write("<fluids>\n")

            file.write(f"\t<Gc> {Gc} </Gc>\n")
            file.write(f"\t<omega_f1> {omega_f1} </omega_f1>\n")
            file.write(f"\t<omega_f2> {omega_f2} </omega_f2>\n")
            file.write(f"\t<force_f1> {force_f1} </force_f1>\n")
            file.write(f"\t<force_f2> {force_f2} </force_f2>\n")
            file.write(f"\t<G_ads_f1_s1> {G_ads_f1_s1} </G_ads_f1_s1>\n")
            file.write(f"\t<G_ads_f1_s2> {G_ads_f1_s2} </G_ads_f1_s2>\n")
            file.write(f"\t<G_ads_f1_s3> {G_ads_f1_s3} </G_ads_f1_s3>\n")
            file.write(f"\t<G_ads_f1_s4> {G_ads_f1_s4} </G_ads_f1_s4>\n")

            file.write(f"\t<rho_f1> {rho_f1} </rho_f1>\n")
            file.write(f"\t<rho_f2> {rho_f2} </rho_f2>\n")

            file.write(f"\t<pressure_bc> {pressure_bc} </pressure_bc>\n")
            file.write(f"\t<rho_f1_i> {rho_f1} </rho_f1_i>\n")
            file.write(f"\t<rho_f2_i> {rho_f2} </rho_f2_i>\n")
            file.write(f"\t<num_pc_steps> {num_pc_steps} </num_pc_steps>\n")
            file.write(f"\t<min_radius> {minimum_radius} </min_radius>\n")
            file.write(f"\t<rho_d> 0.06 </rho_d>\n")

            file.write("</fluids>\n\n")

    def _write_output_section(self):
        convergence = self.inputs["simulation"]["convergence"]
        convergence_iter = self.inputs["simulation"]["convergence_iter"]
        max_iter = self.inputs["simulation"]["max_iterations"]
        save_sim = self.inputs["simulation"]["save_sim"]
        save_iter = self.inputs["simulation"]["save_iter"]
        gif_iter = self.inputs["simulation"]["gif_iter"]
        vtk_iter = self.inputs["simulation"]["vtk_iter"]
        rho_f2_vtk = self.inputs["simulation"]["rho_f2_vtk"]
        print_geom = self.inputs["simulation"]["print_geom"]
        print_stl = self.inputs["simulation"]["print_stl"]

        with open(self.file, "a") as file:
            # Write output section
            file.write("<output>\n")

            # The / after the {dir} is important for the c++ code
            file.write(f"\t<out_folder> {self.output_dir}/ </out_folder>\n")

            file.write(f"\t<save_it> {save_iter} </save_it>\n")
            file.write(f"\t<save_sim> {save_sim} </save_sim>\n")
            file.write(f"\t<convergence> {convergence} </convergence>\n")
            file.write(f"\t<it_max> {max_iter} </it_max>\n")
            file.write(f"\t<it_conv> {convergence_iter} </it_conv>\n")
            file.write(f"\t<it_gif> {gif_iter} </it_gif>\n")
            file.write(f"\t<rho_vtk> {rho_f2_vtk} </rho_vtk>\n")
            file.write(f"\t<it_vtk> {vtk_iter} </it_vtk>\n")
            file.write(f"\t<print_geom> {print_geom} </print_geom>\n")
            file.write(f"\t<print_stl> {print_stl} </print_stl>\n")

            file.write("</output>")

    def create_input_file(self):
        self._write_geometry_section()
        self._write_fluid_positions()
        self._write_fluid_data()
        self._write_output_section()
