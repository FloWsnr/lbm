"""
Creates the xml input file for the palabos simulation.

Creator: Florian Wiesner
Date: 2024-11-08
"""


class PalabosInputFile:
    def __init__(self, config):
        self.config = config
        self.file = (
            self.config["input_output"]["input_folder"] / "2_phase_sim_input.xml"
        )

        restart_sim = config["simulation"]["restart_sim"]
        with open(self.file, "w") as file:
            file.write('<?xml version="1.0"?>\n\n')
            file.write(f"<load_savedstated> {restart_sim} </load_savedstated>\n\n")

    def _write_geometry_section(self):
        # Get geometry inputs
        periodic_x = self.config["domain"]["periodic_boundary"]["x"]
        periodic_y = self.config["domain"]["periodic_boundary"]["y"]
        periodic_z = self.config["domain"]["periodic_boundary"]["z"]

        structure = self.config["structure"]
        nx = structure.shape[0]
        ny = structure.shape[1]
        nz = structure.shape[2]
        domain_size = [nx, ny, nz]

        geo_file_name = self.config["input_output"]["file_name"]
        # Remove suffix
        geo_file_name = geo_file_name.split(".")[0]
        geo_file = self.config["input_output"]["input_folder"] / f"{geo_file_name}.dat"

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
            file.write("\t<per>\n")
            file.write(
                f"\t\t<fluid1> <x> {periodic_x} </x> <y> {periodic_y} </y> <z> {periodic_z} </z> </fluid1>\n"
            )
            file.write(
                f"\t\t<fluid2> <x> {periodic_x} </x> <y> {periodic_y} </y> <z> {periodic_z} </z> </fluid2>\n"
            )
            file.write("\t</per>\n")
            file.write("</geometry>\n\n")

    def _write_fluid_positions(self):
        load_fluid_type = self.config["simulation"]["fluid_init"]
        if load_fluid_type == "geom":
            load_fluid_from_geom = True
        else:
            load_fluid_from_geom = False

        fluid1_x1 = self.config["simulation"]["fluid_1_init"]["x1"]
        fluid1_x2 = self.config["simulation"]["fluid_1_init"]["x2"]
        fluid1_y1 = self.config["simulation"]["fluid_1_init"]["y1"]
        fluid1_y2 = self.config["simulation"]["fluid_1_init"]["y2"]
        fluid1_z1 = self.config["simulation"]["fluid_1_init"]["z1"]
        fluid1_z2 = self.config["simulation"]["fluid_1_init"]["z2"]
        fluid2_x1 = self.config["simulation"]["fluid_2_init"]["x1"]
        fluid2_x2 = self.config["simulation"]["fluid_2_init"]["x2"]
        fluid2_y1 = self.config["simulation"]["fluid_2_init"]["y1"]
        fluid2_y2 = self.config["simulation"]["fluid_2_init"]["y2"]
        fluid2_z1 = self.config["simulation"]["fluid_2_init"]["z1"]
        fluid2_z2 = self.config["simulation"]["fluid_2_init"]["z2"]

        with open(self.file, "a") as file:
            # Write initial position of fluids
            file.write("<init>\n")
            file.write(
                f"\t<fluid_from_geom> {load_fluid_from_geom} </fluid_from_geom>\n"
            )
            file.write("\t<fluid1>\n")
            file.write(
                f"\t\t <x1> {fluid1_x1} </x1> <y1> {fluid1_y1} </y1> <z1> {fluid1_z1} </z1>\n"
            )
            file.write(
                f"\t\t <x2> {fluid1_x2} </x2> <y2> {fluid1_y2} </y2> <z2> {fluid1_z2} </z2>\n"
            )
            file.write("\t</fluid1>\n")
            file.write("\t<fluid2>\n")
            file.write(
                f"\t\t <x1> {fluid2_x1} </x1> <y1> {fluid2_y1} </y1> <z1> {fluid2_z1} </z1>\n"
            )
            file.write(
                f"\t\t <x2> {fluid2_x2} </x2> <y2> {fluid2_y2} </y2> <z2> {fluid2_z2} </z2>\n"
            )
            file.write("\t</fluid2>\n")
            file.write("</init>\n\n")

    def _write_fluid_data(self):
        rho_f1 = self.config["simulation"]["rho_f1"]
        rho_f2 = self.config["simulation"]["rho_f2"]

        force_f1 = self.config["simulation"]["force_f1"]
        force_f2 = self.config["simulation"]["force_f2"]

        Gc = self.config["simulation"]["fluid_data"]["Gc"]
        omega_f1 = self.config["simulation"]["fluid_data"]["omega_f1"]
        omega_f2 = self.config["simulation"]["fluid_data"]["omega_f2"]
        G_ads_f1_s1 = self.config["simulation"]["fluid_data"]["G_ads_f1_s1"]
        G_ads_f1_s2 = self.config["simulation"]["fluid_data"]["G_ads_f1_s2"]
        G_ads_f1_s3 = self.config["simulation"]["fluid_data"]["G_ads_f1_s3"]
        G_ads_f1_s4 = self.config["simulation"]["fluid_data"]["G_ads_f1_s4"]

        pressure_bc = self.config["simulation"]["pressure_bc"]
        if pressure_bc:
            minimum_radius = self.config["simulation"]["minimum_radius"]
            num_pc_steps = self.config["simulation"]["num_pressure_steps"]
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
            file.write("\t<rho_d> 0.06 </rho_d>\n")

            file.write("</fluids>\n\n")

    def _write_output_section(self):
        convergence = self.config["simulation"]["convergence"]
        convergence_iter = self.config["simulation"]["convergence_iter"]
        max_iter = self.config["simulation"]["max_iter"]
        max_iter_per_p = self.config["simulation"]["max_iter_per_p"]
        save_sim = self.config["simulation"]["save_sim"]
        save_iter = self.config["simulation"]["save_iter"]
        gif_iter = self.config["simulation"]["gif_iter"]
        vtk_iter = self.config["simulation"]["vtk_iter"]
        rho_f2_vtk = self.config["simulation"]["rho_f2_vtk"]
        print_geom = self.config["simulation"]["print_geom"]
        print_stl = self.config["simulation"]["print_stl"]

        with open(self.file, "a") as file:
            # Write output section
            file.write("<output>\n")

            # The / after the {dir} is important for the c++ code
            output_dir = self.config["input_output"]["output_folder"]
            file.write(f"\t<out_folder> {output_dir}/ </out_folder>\n")

            file.write(f"\t<save_it> {save_iter} </save_it>\n")
            file.write(f"\t<save_sim> {save_sim} </save_sim>\n")
            file.write(f"\t<convergence> {convergence} </convergence>\n")
            file.write(f"\t<max_iter> {max_iter} </max_iter>\n")
            file.write(f"\t<it_max_per_p> {max_iter_per_p} </it_max_per_p>\n")

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
