"""
Structural analyzer with automatic Core detection.
"""

from pathlib import Path
from typing import List, Optional, Any, Tuple

import numpy as np
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.Structure import Structure
from sklearn.decomposition import PCA

from ..models import CoreResidueInfo, RMSDResult


class StructureAnalyzer:
    """
    Structural analyzer with automatic Core detection.

    Implements intelligent RMSD analysis considering:
    - pLDDT as proxy for structural order
    - Separation between structured Core and flexible IDRs
    - Optimized alignment for fold preservation
    """

    def __init__(self, plddt_threshold: float = 70.0):
        """
        Initialize the analyzer.

        Args:
            plddt_threshold: Minimum pLDDT threshold to consider
                             a residue as part of the Core (default: 70).
        """
        self.plddt_threshold = plddt_threshold
        self.parser = PDBParser(QUIET=True)
        self.superimposer = Superimposer()
        self._wt_structure: Optional[Structure] = None
        self._core_info: Optional[CoreResidueInfo] = None

    def load_wt_structure(self, pdb_path: Path) -> Structure:
        """
        Load and store the Wild Type structure.

        Args:
            pdb_path: Path to WT PDB file (AlphaFold).

        Returns:
            BioPython Structure.
        """
        self._wt_structure = self.parser.get_structure("WT", str(pdb_path))
        print(f"✅ WT structure loaded: {pdb_path.name}")
        return self._wt_structure

    def detect_core_residues(self, pdb_path: Optional[Path] = None) -> CoreResidueInfo:
        """
        Detect structured Core residues based on pLDDT.

        In AlphaFold, pLDDT is stored in the B-factor column of the PDB.

        Args:
            pdb_path: Path to PDB. If None, uses loaded WT structure.

        Returns:
            CoreResidueInfo with indices and values of Core residues.

        Raises:
            ValueError: If no structure is available.
        """
        if pdb_path is not None:
            structure = self.parser.get_structure("temp", str(pdb_path))
        elif self._wt_structure is not None:
            structure = self._wt_structure
        else:
            raise ValueError("No structure loaded. Use load_wt_structure() first.")

        core_indices: List[int] = []
        plddt_values: List[float] = []
        all_plddt: List[float] = []

        # Iterate over residues
        for model in structure:
            for chain in model:
                for res_idx, residue in enumerate(chain.get_residues()):
                    # Skip heteroatoms (water, ligands)
                    if residue.id[0] != ' ':
                        continue

                    # Get pLDDT from CA (stored as B-factor)
                    if 'CA' in residue:
                        plddt = residue['CA'].get_bfactor()
                        all_plddt.append(plddt)

                        if plddt >= self.plddt_threshold:
                            core_indices.append(res_idx)
                            plddt_values.append(plddt)
                break  # First chain only
            break  # First model only

        self._core_info = CoreResidueInfo(
            residue_indices=core_indices,
            plddt_values=plddt_values,
            threshold=self.plddt_threshold
        )

        total_residues = len(all_plddt)
        core_fraction = len(core_indices) / total_residues * 100 if total_residues > 0 else 0

        print(f"\n📊 Auto-Core Detection (pLDDT ≥ {self.plddt_threshold}):")
        print(f"   Total residues: {total_residues}")
        print(f"   Core residues: {len(core_indices)} ({core_fraction:.1f}%)")
        print(f"   Mean pLDDT (Core): {np.mean(plddt_values):.1f}")
        print(f"   Mean pLDDT (Global): {np.mean(all_plddt):.1f}")

        return self._core_info

    def _get_ca_atoms(
        self,
        structure: Structure,
        indices: Optional[List[int]] = None
    ) -> List:
        """
        Extract CA atoms from a structure.

        Args:
            structure: BioPython Structure.
            indices: List of residue indices to extract. If None, gets all.

        Returns:
            List of CA atoms.
        """
        ca_atoms = []

        for model in structure:
            for chain in model:
                residues = list(res for res in chain.get_residues() if res.id[0] == ' ')

                if indices is None:
                    target_residues = residues
                else:
                    target_residues = [residues[i] for i in indices if i < len(residues)]

                for residue in target_residues:
                    if 'CA' in residue:
                        ca_atoms.append(residue['CA'])
                break
            break

        return ca_atoms

    def calculate_rmsd(
        self,
        sample_path: Path,
        sample_id: Optional[str] = None
    ) -> RMSDResult:
        """
        Calculate Core and Global RMSD for a conformation.

        Procedure:
        1. Load sample structure
        2. Extract Core CAs (pre-detected indices)
        3. Align sample to WT using Core only
        4. Calculate Core RMSD
        5. Calculate Global RMSD (with Core transformation applied)

        Args:
            sample_path: Path to conformation PDB.
            sample_id: Sample identifier (optional).

        Returns:
            RMSDResult with calculated metrics.
        """
        if self._wt_structure is None:
            raise ValueError("WT structure not loaded. Use load_wt_structure() first.")
        if self._core_info is None:
            raise ValueError("Core not detected. Use detect_core_residues() first.")

        sample_id = sample_id or sample_path.stem

        # Load sample structure
        sample_structure = self.parser.get_structure(sample_id, str(sample_path))

        # Extract Core CAs
        wt_core_cas = self._get_ca_atoms(self._wt_structure, self._core_info.residue_indices)
        sample_core_cas = self._get_ca_atoms(sample_structure, self._core_info.residue_indices)

        # Check compatibility
        min_core_len = min(len(wt_core_cas), len(sample_core_cas))
        if min_core_len == 0:
            raise ValueError(f"No Core CA atoms found in {sample_id}.")

        wt_core_cas = wt_core_cas[:min_core_len]
        sample_core_cas = sample_core_cas[:min_core_len]

        # Align by Core
        self.superimposer.set_atoms(wt_core_cas, sample_core_cas)
        rmsd_core = self.superimposer.rms

        # Apply transformation to ENTIRE sample structure
        self.superimposer.apply(sample_structure.get_atoms())

        # Calculate Global RMSD (after Core alignment)
        wt_all_cas = self._get_ca_atoms(self._wt_structure)
        sample_all_cas = self._get_ca_atoms(sample_structure)

        min_total_len = min(len(wt_all_cas), len(sample_all_cas))
        wt_all_cas = wt_all_cas[:min_total_len]
        sample_all_cas = sample_all_cas[:min_total_len]

        # Calculate Global RMSD manually (structure already aligned)
        wt_coords = np.array([atom.get_coord() for atom in wt_all_cas])
        sample_coords = np.array([atom.get_coord() for atom in sample_all_cas])

        diff = wt_coords - sample_coords
        rmsd_global = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

        return RMSDResult(
            sample_id=sample_id,
            sample_path=sample_path,
            rmsd_core=rmsd_core,
            rmsd_global=rmsd_global,
            num_core_atoms=min_core_len,
            num_total_atoms=min_total_len
        )

    def calculate_pca(
        self,
        sample_paths: List[Path]
    ) -> Tuple[np.ndarray, List[str], PCA]:
        """
        Analyse the conformations and generate a PCA.

        Args:
            sample_paths: List of paths to conformation PDBs.

        Returns:
            Tuple of (principal_components, sample_ids, pca_model)
        """
        pca_data = []
        sample_ids = []
        
        for i, path in enumerate(sample_paths):
            try:
                sample_id = f"conf_{i:03d}"

                # Get the structure
                sample_structure = self.parser.get_structure(sample_id, str(path))

                # Extract Core CAs
                wt_core_cas = self._get_ca_atoms(self._wt_structure, self._core_info.residue_indices)
                sample_core_cas = self._get_ca_atoms(sample_structure, self._core_info.residue_indices)

                # Check compatibility
                min_core_len = min(len(wt_core_cas), len(sample_core_cas))
                if min_core_len == 0:
                    raise ValueError(f"No Core CA atoms found in {sample_id}.")

                wt_core_cas = wt_core_cas[:min_core_len]
                sample_core_cas = sample_core_cas[:min_core_len]

                # Align by Core
                self.superimposer.set_atoms(wt_core_cas, sample_core_cas)

                # Apply transformation to ENTIRE sample structure
                self.superimposer.apply(sample_structure.get_atoms())

                wt_all_cas = self._get_ca_atoms(self._wt_structure)
                sample_all_cas = self._get_ca_atoms(sample_structure)

                min_total_len = min(len(wt_all_cas), len(sample_all_cas))
                sample_all_cas = sample_all_cas[:min_total_len]

                sample_coords = np.array([atom.get_coord() for atom in sample_all_cas])
                flatten_sample = sample_coords.flatten()

                pca_data.append(flatten_sample)
                sample_ids.append(sample_id)

            except Exception as e:
                print(f"   ⚠️  Error analyzing {path.name}: {e}")

        X = np.array(pca_data)

        n_samples = X.shape[0]
        if n_samples < 2:
            raise ValueError("Insufficient samples for PCA.")

        pca_model = PCA(n_components=2)
        principal_components = pca_model.fit_transform(X)

        return principal_components, sample_ids, pca_model

    def analyze_all_samples(
        self, 
        wt_paths: List[Path], 
        sample_paths: List[Path]
    ) -> List[RMSDResult]:
        """
        Analyze all generated conformations.

        Args:
            wt_paths: List of paths to WT PDBs.
            sample_paths: List of paths to conformation PDBs.

        Returns:
            List of RMSDResults sorted by Core RMSD.
        """
        print(f"\n🔍 Analyzing {len(sample_paths)} conformations...")

        results: List[RMSDResult] = []

        for i, path in enumerate(sample_paths):
            try:
                result = self.calculate_rmsd(path, sample_id=f"conf_{i:03d}")
                results.append(result)
                print(f"   📐 {result.sample_id}: Core={result.rmsd_core:.3f}Å, Global={result.rmsd_global:.3f}Å")
            except Exception as e:
                print(f"   ⚠️  Error analyzing {path.name}: {e}")

        # Sort by Core RMSD
        results.sort(key=lambda x: x.rmsd_core)

        return results

    def get_best_candidate(self, results: List[RMSDResult]) -> RMSDResult:
        """
        Select the best conformation for RFdiffusion.

        Criterion: Lowest Core RMSD (maximum structural fold preservation).

        Args:
            results: List of RMSDResults.

        Returns:
            RMSDResult of the best conformation.
        """
        if not results:
            raise ValueError("Empty results list.")

        best = min(results, key=lambda x: x.rmsd_core)

        print(f"\n🏆 Best candidate for RFdiffusion:")
        print(f"   ID: {best.sample_id}")
        print(f"   Core RMSD: {best.rmsd_core:.3f} Å")
        print(f"   Global RMSD: {best.rmsd_global:.3f} Å")
        print(f"   Flexibility Score: {best.flexibility_score:.2f}")

        return best

    @property
    def core_info(self) -> Optional[CoreResidueInfo]:
        """Get the detected core residue information."""
        return self._core_info

    @property
    def wt_structure(self) -> Optional[Structure]:
        """Get the loaded WT structure."""
        return self._wt_structure
