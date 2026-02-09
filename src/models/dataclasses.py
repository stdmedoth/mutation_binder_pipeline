"""
Dataclasses for mutation information and analysis results.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List

import numpy as np

from ..constants import AA_CODE_MAP


@dataclass
class MutationInfo:
    """Data structure for mutation information."""
    original_aa: str  # Original amino acid (1-letter)
    position: int     # Position in sequence (1-indexed)
    mutant_aa: str    # Mutant amino acid (1-letter)

    @property
    def notation(self) -> str:
        """Returns standard mutation notation (e.g., D40G)."""
        return f"{self.original_aa}{self.position}{self.mutant_aa}"

    @classmethod
    def from_string(cls, mutation_str: str) -> 'MutationInfo':
        """
        Parses mutation string (e.g., 'D40G').
        
        Args:
            mutation_str: Mutation in format [original_aa][position][mutant_aa]
            
        Returns:
            MutationInfo instance
            
        Raises:
            ValueError: If mutation format is invalid
        """
        mutation_str = mutation_str.strip().upper()

        if len(mutation_str) < 3:
            raise ValueError(f"Invalid mutation format: '{mutation_str}'. Use format like 'D40G'.")

        original_aa = mutation_str[0]
        mutant_aa = mutation_str[-1]

        try:
            position = int(mutation_str[1:-1])
        except ValueError:
            raise ValueError(f"Invalid position in mutation: '{mutation_str}'.")

        if original_aa not in AA_CODE_MAP:
            raise ValueError(f"Invalid original amino acid: '{original_aa}'.")
        if mutant_aa not in AA_CODE_MAP:
            raise ValueError(f"Invalid mutant amino acid: '{mutant_aa}'.")
        if position < 1:
            raise ValueError(f"Position must be >= 1, received: {position}.")

        return cls(original_aa=original_aa, position=position, mutant_aa=mutant_aa)


@dataclass
class CoreResidueInfo:
    """Information about structured Core residues."""
    residue_indices: List[int]  # 0-based indices of core residues
    plddt_values: List[float]   # Corresponding pLDDT values
    threshold: float            # Threshold used for filtering

    @property
    def num_residues(self) -> int:
        """Number of core residues."""
        return len(self.residue_indices)

    @property
    def mean_plddt(self) -> float:
        """Mean pLDDT of core residues."""
        return np.mean(self.plddt_values) if self.plddt_values else 0.0


@dataclass
class RMSDResult:
    """RMSD analysis result."""
    sample_id: str
    sample_path: Path
    rmsd_core: float       # RMSD of Core residues only
    rmsd_global: float     # RMSD of all residues (after Core alignment)
    num_core_atoms: int    # Number of CA atoms in Core
    num_total_atoms: int   # Total number of CA atoms

    @property
    def flexibility_score(self) -> float:
        """Flexibility score: RMSD_global/RMSD_core ratio."""
        if self.rmsd_core > 0:
            return self.rmsd_global / self.rmsd_core
        return float('inf')
