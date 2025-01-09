"""
Post processing of the vti files to extract saturation and pressure.
Store the results in a csv file.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Union

from lbm_wetting.utils.postprocessing.vti_processing import read_vti_file


class SaturationCalculator:
    def __init__(self, processed_dir: Union[Path, str]) -> None:
        """Initialize the SaturationCalculator.

        Parameters
        ----------
        processed_dir : Union[Path, str]
            Directory containing the processed VTI files
        """
        self.processed_dir = Path(processed_dir)
        if not self.processed_dir.exists():
            raise FileNotFoundError(f"Directory {processed_dir} does not exist")

    def _find_vti_files(self) -> list[Path]:
        """Find all VTI files in the processed directory."""
        vti_files = list(self.processed_dir.glob("sim_*.vti"))
        # Sort by simulation step and time
        vti_files.sort(
            key=lambda x: (int(x.stem.split("_")[1]), int(x.stem.split("_")[2]))
        )
        return vti_files

    def _calculate_saturation(self, vti_file: Path) -> float:
        """Calculate saturation (percentage of non-zero pixels) for a single VTI file.

        Parameters
        ----------
        vti_file : Path
            Path to the VTI file

        Returns
        -------
        float
            Saturation value (between 0 and 1)
        """
        data = read_vti_file(vti_file)
        wet_structure = data["wet_structure"]

        # Count water pixels
        total_pixels = wet_structure.size
        water_pixels = np.count_nonzero(wet_structure == 128)

        return water_pixels / total_pixels

    def _load_pressure_diffs(self, output_file: Path) -> pd.DataFrame:
        """Load the pressure differences from the output file."""
        with open(output_file, "r") as file:
            lines = file.readlines()

        # Extract pressure difference lines
        pressure_lines = [line for line in lines if "Pressure difference" in line]
        if not pressure_lines:
            raise ValueError("No pressure difference lines found in the output file")

        # Extract pressure differences
        pressure_diffs = [float(line.split()[-1]) for line in pressure_lines]

        return pressure_diffs

    def process_all_files(self) -> pd.DataFrame:
        """Process all VTI files and calculate saturations.

        Returns
        -------
        pd.DataFrame
            DataFrame containing timesteps and corresponding saturation values
        """
        vti_files = self._find_vti_files()
        output_file = self.processed_dir.parent / "output" / "output.dat"
        pressure_diffs = self._load_pressure_diffs(output_file)
        results = []

        for file in vti_files:
            # Extract timestep from filename (sim_001_00002000.vti -> 1.2000)
            step = int(file.stem.split("_")[1])
            iter = int(file.stem.split("_")[2])

            saturation = self._calculate_saturation(file)
            results.append(
                {
                    "pressure": pressure_diffs[step],
                    "step": step,
                    "iteration": iter,
                    "saturation": saturation,
                }
            )

        df = pd.DataFrame(results)
        return df

    def save_results(self, df: pd.DataFrame) -> None:
        """Calculate saturations and save results to CSV.

        Parameters
        ----------
        output_file : Union[Path, str], optional
            Path to save the CSV file. If None, saves to 'saturation.csv' in processed directory
        """
        output_file = self.processed_dir / "saturation.csv"
        df.to_csv(output_file, index=False)
        print(f"Saved saturation results to {output_file}")


if __name__ == "__main__":
    import argparse

    sim_dir = Path("/hpcwork/fw641779/lbm/Toray-120C/55cov/structure0")
    sim_name = "test_run_2"

    parser = argparse.ArgumentParser(description="Calculate saturation from VTI files")
    parser.add_argument(
        "--sim-dir",
        "-s",
        type=str,
        help="Directory containing simulation data",
        default=sim_dir,
    )
    parser.add_argument(
        "--sim-name",
        "-n",
        type=str,
        help="Name of the simulation",
        default=sim_name,
    )
    args = parser.parse_args()

    processed_dir = Path(args.sim_dir) / args.sim_name / "processed"

    calculator = SaturationCalculator(processed_dir)
    df = calculator.process_all_files()
    calculator.save_results(df)
