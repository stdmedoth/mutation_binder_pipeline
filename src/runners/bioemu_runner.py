"""
Wrapper for BioEmu conformational sampling model.
"""

from pathlib import Path
from typing import List, Optional, Callable

import numpy as np

from ..constants import AA_CODE_MAP


class BioEmuRunner:
    """
    Wrapper for the BioEmu conformational sampling model.

    BioEmu uses a diffusion model to generate protein conformations
    that sample the Boltzmann ensemble of the protein.
    """

    def __init__(self, output_dir: Path):
        """
        Initialize the BioEmu runner.

        Args:
            output_dir: Directory to save generated conformations.
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.sample_func: Optional[Callable] = None
        self._initialized = False

    def initialize_model(self) -> None:
        """
        Load the BioEmu model from the correct submodule.
        """
        if self._initialized:
            print("📋 BioEmu model already initialized.")
            return

        try:
            import torch
            if not torch.cuda.is_available():
                raise RuntimeError("GPU not available. BioEmu requires CUDA.")

            print("🔄 Loading BioEmu model...")

            # CORRECT IMPORT: The sampling function is in bioemu.sample
            from bioemu.sample import main as bioemu_sample

            self.sample_func = bioemu_sample
            self._initialized = True
            print(f"✅ BioEmu initialized successfully (GPU: {torch.cuda.get_device_name(0)})")

        except ImportError as e:
            raise RuntimeError(
                f"BioEmu import error: {e}\n"
                "Ensure you have installed the package: pip install bioemu"
            ) from e

    def run_sampling(
        self,
        sequence: str,
        num_samples: int = 10,
        prefix: str = "mutant"
    ) -> List[Path]:
        """
        Generate structural conformations for a sequence.
        
        Args:
            sequence: Amino acid sequence.
            num_samples: Number of conformations to generate.
            prefix: Prefix for output directory.
            
        Returns:
            List of paths to generated PDB files.
        """
        if not self._initialized:
            self.initialize_model()

        print(f"\n🔬 Generating {num_samples} conformations for {len(sequence)} aa sequence...")

        # Create a specific subdirectory for this run to avoid file conflicts
        run_dir = self.output_dir / prefix
        run_dir.mkdir(parents=True, exist_ok=True)

        try:
            # Run sampling
            self.sample_func(
                sequence=sequence,
                num_samples=num_samples,
                output_dir=str(run_dir)
            )

            # Collect the generated files
            generated_files = sorted(list(run_dir.glob("*.pdb")))

            if not generated_files:
                raise RuntimeError("BioEmu finished but no PDB files were found in output directory.")

            print(f"\n✅ {len(generated_files)} conformations generated in: {run_dir}")
            return generated_files

        except Exception as e:
            raise RuntimeError(f"Error during BioEmu sampling: {e}") from e

    def _coords_to_pdb(
        self,
        coords: np.ndarray,
        sequence: str,
        output_path: Path
    ) -> None:
        """
        Convert backbone coordinates to PDB format.

        Args:
            coords: Coordinate array. Expected shapes:
                    - (N, 3) for CA-only
                    - (N, 4, 3) for N, CA, C, O
                    - (N, 37, 3) for all-atom
            sequence: Amino acid sequence.
            output_path: Output path for PDB.
        """
        coords = np.array(coords)

        # Determine coordinate type based on shape
        if coords.ndim == 2:
            # (N, 3) - CA only
            atom_names = ['CA']
            coords = coords[:, np.newaxis, :]  # Make (N, 1, 3)
        elif coords.ndim == 3:
            if coords.shape[1] == 4:
                atom_names = ['N', 'CA', 'C', 'O']
            elif coords.shape[1] == 3:
                atom_names = ['N', 'CA', 'C']
            elif coords.shape[1] == 37:
                # All-atom representation, extract backbone
                atom_names = ['N', 'CA', 'C', 'O']
                coords = coords[:, :4, :]
            else:
                atom_names = ['CA']
                coords = coords[:, :1, :]
        else:
            raise ValueError(f"Unexpected coordinate shape: {coords.shape}")

        with open(output_path, 'w') as f:
            atom_serial = 1

            for res_idx, aa in enumerate(sequence):
                res_name = AA_CODE_MAP.get(aa, 'UNK')
                res_num = res_idx + 1

                for atom_idx, atom_name in enumerate(atom_names):
                    if res_idx < len(coords) and atom_idx < coords.shape[1]:
                        x, y, z = coords[res_idx, atom_idx]

                        # Skip if coordinates are NaN or zero
                        if np.isnan(x) or np.isnan(y) or np.isnan(z):
                            continue

                        f.write(
                            f"ATOM  {atom_serial:5d}  {atom_name:<3s} {res_name:3s} A"
                            f"{res_num:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
                            f"  1.00  0.00           {atom_name[0]:>2s}\n"
                        )
                        atom_serial += 1

            f.write("END\n")
