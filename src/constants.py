"""
Constants and amino acid mappings for the mutation pipeline.
"""

from typing import Dict

# Amino acid code mapping (1-letter <-> 3-letter)
AA_CODE_MAP: Dict[str, str] = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
}

# Reverse mapping (3-letter -> 1-letter)
AA_3TO1: Dict[str, str] = {v: k for k, v in AA_CODE_MAP.items()}

# API URLs
UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
ALPHAFOLD_PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v6.pdb"
