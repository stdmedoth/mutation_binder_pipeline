"""
Structure alignment utilities.
"""

from pathlib import Path
from typing import Tuple

from Bio.PDB import PDBParser, PDBIO, Superimposer


def align_and_save_best_candidate(pipeline) -> Path:
    """
    Align the best candidate to WT (using Core) and save to disk.
    Necessary for correct visualization and for RFdiffusion.
    
    Args:
        pipeline: MutationPipeline instance with results
        
    Returns:
        Path to aligned PDB file
    """
    print("\n🔄 Aligning best candidate for visualization...")

    # 1. Load structures
    parser = PDBParser(QUIET=True)
    wt_struct = parser.get_structure("WT", str(pipeline.wt_structure_path))
    mutant_struct = parser.get_structure("Mutant", str(pipeline.best_result.sample_path))

    # 2. Get Core atoms (based on previously detected indices)
    core_indices = pipeline.analyzer.core_info.residue_indices

    def get_core_atoms(structure, indices):
        atoms = []
        for model in structure:
            for chain in model:
                residues = list(chain.get_residues())
                for i in indices:
                    if i < len(residues) and 'CA' in residues[i]:
                        atoms.append(residues[i]['CA'])
                return atoms  # Return after first chain

    wt_core = get_core_atoms(wt_struct, core_indices)
    mutant_core = get_core_atoms(mutant_struct, core_indices)

    # 3. Calculate and apply rotation matrix
    si = Superimposer()
    min_len = min(len(wt_core), len(mutant_core))
    si.set_atoms(wt_core[:min_len], mutant_core[:min_len])
    si.apply(mutant_struct.get_atoms())  # Apply rotation to ALL atoms of mutant

    print(f"   Alignment RMSD: {si.rms:.3f} Å")

    # 4. Save aligned file
    output_path = pipeline.best_result.sample_path.parent / f"{pipeline.best_result.sample_id}_aligned.pdb"

    io = PDBIO()
    io.set_structure(mutant_struct)
    io.save(str(output_path))

    print(f"✅ Aligned file saved: {output_path.name}")
    
    return output_path


def align_structures(
    reference_path: Path,
    mobile_path: Path,
    core_indices: list,
    output_path: Path = None
) -> Tuple[Path, float]:
    """
    Align mobile structure to reference using Core residues.
    
    Args:
        reference_path: Path to reference PDB
        mobile_path: Path to mobile PDB
        core_indices: List of residue indices for alignment
        output_path: Output path for aligned structure (optional)
        
    Returns:
        Tuple of (output_path, rmsd)
    """
    parser = PDBParser(QUIET=True)
    ref_struct = parser.get_structure("ref", str(reference_path))
    mobile_struct = parser.get_structure("mobile", str(mobile_path))

    def get_ca_atoms(structure, indices):
        atoms = []
        for model in structure:
            for chain in model:
                residues = list(chain.get_residues())
                for i in indices:
                    if i < len(residues) and 'CA' in residues[i]:
                        atoms.append(residues[i]['CA'])
                return atoms

    ref_atoms = get_ca_atoms(ref_struct, core_indices)
    mobile_atoms = get_ca_atoms(mobile_struct, core_indices)

    si = Superimposer()
    min_len = min(len(ref_atoms), len(mobile_atoms))
    si.set_atoms(ref_atoms[:min_len], mobile_atoms[:min_len])
    si.apply(mobile_struct.get_atoms())

    if output_path is None:
        output_path = mobile_path.parent / f"{mobile_path.stem}_aligned.pdb"

    io = PDBIO()
    io.set_structure(mobile_struct)
    io.save(str(output_path))

    return output_path, si.rms
