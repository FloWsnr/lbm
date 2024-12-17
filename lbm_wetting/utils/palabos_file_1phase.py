"""
Creates the xml input file for the palabos simulation.

Creator: Florian Wiesner
Date: 2024-11-08
"""


class PalabosInputFile1Phase:
    def __init__(self, config):
        self.config = config
        self.file = (
            self.config["input_output"]["input_folder"] / "1_phase_sim_input.xml"
        )

        with open(self.file, "w") as file:
            file.write('<?xml version="1.0"?>\n\n')

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

    def _write_output_section(self):
        convergence = self.config["simulation"]["convergence"]
        max_iter = self.config["simulation"]["max_iter"]
        pressure = self.config["simulation"]["pressure"]
        output_dir = self.config["input_output"]["output_folder"]
        input_dir = self.config["input_output"]["input_folder"]

        with open(self.file, "a") as file:
            # Write output section
            file.write("<folder>\n")

            # The / after the {dir} is important for the c++ code
            file.write(f"\t<out_f> {output_dir}/ </out_f>\n")
            file.write(f"\t<in_f> {input_dir}/ </in_f>\n")

            file.write("</folder>")

            file.write("</simulations>\n")
            file.write(f"\t<press> {pressure} </press>\n")
            file.write("\t<num> 1 </num>\n")
            file.write(f"\t<iter> {max_iter} </iter>\n")
            file.write(f"\t<conv> {convergence} </conv>\n")
            file.write("\t<vtk_out> True </vtk_out>\n")
            file.write("</simulations>")

    def create_input_file(self):
        self._write_geometry_section()
        self._write_output_section()
