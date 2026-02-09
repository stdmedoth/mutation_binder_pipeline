"""
Clean PDB exporter for RFdiffusion.
"""

from pathlib import Path
from typing import List

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select

from ..models import MutationInfo, CoreResidueInfo, RMSDResult


class MutantPDBExporter:
    """
    Clean PDB exporter for RFdiffusion.
    """

    @staticmethod
    def export_clean_pdb(
        source_path: Path,
        output_path: Path,
        remove_hydrogens: bool = True,
        remove_water: bool = True
    ) -> Path:
        """
        Export clean PDB ready for RFdiffusion.

        Args:
            source_path: Source PDB path.
            output_path: Output path.
            remove_hydrogens: Remove hydrogen atoms.
            remove_water: Remove water molecules.

        Returns:
            Path of exported file.
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("mutant", str(source_path))

        class CleanSelect(Select):
            def accept_atom(self, atom):
                # Remove hydrogens
                if remove_hydrogens and atom.element == 'H':
                    return False
                return True

            def accept_residue(self, residue):
                # Remove water
                if remove_water and residue.id[0] == 'W':
                    return False
                # Accept only standard residues
                return residue.id[0] == ' '

        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_path), CleanSelect())

        print(f"✅ Clean PDB exported: {output_path}")
        return output_path

    @staticmethod
    def generate_summary(
        uniprot_id: str,
        mutation: MutationInfo,
        best_result: RMSDResult,
        all_results: List[RMSDResult],
        core_info: CoreResidueInfo
    ) -> str:
        """
        Generate executive summary of the analysis.

        Returns:
            Formatted string with summary.
        """
        summary = f"""
╔══════════════════════════════════════════════════════════════════════╗
║              EXECUTIVE SUMMARY - MUTATION PIPELINE                   ║
╠══════════════════════════════════════════════════════════════════════╣
║                                                                      ║
║  📋 PROTEIN INFORMATION                                              ║
║  ├─ UniProt ID: {uniprot_id:<52} ║
║  └─ Mutation: {mutation.notation:<55} ║
║                                                                      ║
║  🔬 STRUCTURAL ANALYSIS                                              ║
║  ├─ pLDDT Threshold: {core_info.threshold:<47.0f} ║
║  ├─ Core Residues: {core_info.num_residues:<49} ║
║  └─ Mean Core pLDDT: {core_info.mean_plddt:<46.1f} ║
║                                                                      ║
║  📊 RMSD STATISTICS ({len(all_results)} conformations)                             ║
║  ├─ Core RMSD (mean): {np.mean([r.rmsd_core for r in all_results]):<45.3f} ║
║  ├─ Core RMSD (std): {np.std([r.rmsd_core for r in all_results]):<46.3f} ║
║  ├─ Global RMSD (mean): {np.mean([r.rmsd_global for r in all_results]):<42.3f} ║
║  └─ Global RMSD (std): {np.std([r.rmsd_global for r in all_results]):<44.3f} ║
║                                                                      ║
║  🏆 BEST CANDIDATE FOR RFdiffusion                                   ║
║  ├─ ID: {best_result.sample_id:<59} ║
║  ├─ Core RMSD: {best_result.rmsd_core:<53.3f} ║
║  ├─ Global RMSD: {best_result.rmsd_global:<51.3f} ║
║  └─ Flexibility Score: {best_result.flexibility_score:<44.2f} ║
║                                                                      ║
╚══════════════════════════════════════════════════════════════════════╝
"""
        return summary
