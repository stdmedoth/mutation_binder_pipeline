"""
Handler for protein data retrieval and manipulation.
"""

from pathlib import Path
from typing import Dict, Tuple

import requests

from ..constants import AA_CODE_MAP, UNIPROT_FASTA_URL, ALPHAFOLD_PDB_URL
from ..models import MutationInfo


class ProteinDataHandler:
    """
    Handler for protein data retrieval and manipulation.

    Responsibilities:
    - Fetch FASTA sequences from UniProt
    - Download AlphaFold structures
    - Apply and validate mutations
    """

    def __init__(self, output_dir: Path):
        """
        Initialize the handler.

        Args:
            output_dir: Directory to save downloaded structures.
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._sequence_cache: Dict[str, str] = {}

    def fetch_sequence(self, uniprot_id: str) -> str:
        """
        Fetch official FASTA sequence from UniProt.

        Args:
            uniprot_id: UniProt identifier (e.g., 'P12345').

        Returns:
            Amino acid sequence as string.

        Raises:
            ValueError: If ID not found or request error.
        """
        uniprot_id = uniprot_id.strip().upper()

        # Check cache
        if uniprot_id in self._sequence_cache:
            print(f"📋 Sequence {uniprot_id} retrieved from cache.")
            return self._sequence_cache[uniprot_id]

        url = UNIPROT_FASTA_URL.format(uniprot_id=uniprot_id)

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            if response.status_code == 404:
                raise ValueError(f"UniProt ID '{uniprot_id}' not found.") from e
            raise ValueError(f"Error fetching sequence: {e}") from e
        except requests.exceptions.RequestException as e:
            raise ValueError(f"Connection error to UniProt: {e}") from e

        # Parse FASTA
        lines = response.text.strip().split('\n')
        sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))

        if not sequence:
            raise ValueError(f"Empty sequence returned for {uniprot_id}.")

        # Validate sequence
        invalid_chars = set(sequence) - set(AA_CODE_MAP.keys())
        if invalid_chars:
            print(f"⚠️  Non-standard characters found and will be ignored: {invalid_chars}")
            sequence = ''.join(aa for aa in sequence if aa in AA_CODE_MAP)

        self._sequence_cache[uniprot_id] = sequence
        print(f"✅ Sequence {uniprot_id} downloaded: {len(sequence)} amino acids.")

        return sequence

    def fetch_alphafold_structure(self, uniprot_id: str) -> Path:
        """
        Download predicted structure from AlphaFold Database (v6).

        Args:
            uniprot_id: UniProt identifier.

        Returns:
            Path to downloaded PDB file.

        Raises:
            ValueError: If structure not available.
        """
        uniprot_id = uniprot_id.strip().upper()
        output_path = self.output_dir / f"AF-{uniprot_id}-F1-model_v6.pdb"

        # Check if already exists
        if output_path.exists():
            print(f"📂 AlphaFold structure already exists: {output_path}")
            return output_path

        url = ALPHAFOLD_PDB_URL.format(uniprot_id=uniprot_id)

        try:
            response = requests.get(url, timeout=60)
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            if response.status_code == 404:
                raise ValueError(
                    f"AlphaFold structure not available for '{uniprot_id}'. "
                    f"Please verify the UniProt ID is correct."
                ) from e
            raise ValueError(f"Error downloading structure: {e}") from e
        except requests.exceptions.RequestException as e:
            raise ValueError(f"Connection error to AlphaFold DB: {e}") from e

        # Save file
        with open(output_path, 'w') as f:
            f.write(response.text)

        print(f"✅ AlphaFold structure downloaded: {output_path}")
        return output_path

    def apply_mutation(
        self,
        sequence: str,
        mutation_str: str
    ) -> Tuple[str, MutationInfo]:
        """
        Apply a point mutation to the sequence.

        Args:
            sequence: Original amino acid sequence.
            mutation_str: Mutation string in format 'D40G'
                          (Asp at position 40 -> Gly).

        Returns:
            Tuple containing:
            - Mutated sequence
            - MutationInfo object with mutation details

        Raises:
            ValueError: If mutation is invalid or inconsistent.
        """
        mutation = MutationInfo.from_string(mutation_str)

        # Validate position
        if mutation.position > len(sequence):
            raise ValueError(
                f"Position {mutation.position} exceeds sequence length ({len(sequence)} aa)."
            )

        # Validate original amino acid (position 1-indexed -> 0-indexed)
        idx = mutation.position - 1
        actual_aa = sequence[idx]

        if actual_aa != mutation.original_aa:
            raise ValueError(
                f"Mutation inconsistency {mutation.notation}: "
                f"expected '{mutation.original_aa}' at position {mutation.position}, "
                f"but found '{actual_aa}'."
            )

        # Apply mutation
        mutated_sequence = sequence[:idx] + mutation.mutant_aa + sequence[idx + 1:]

        print(f"✅ Mutation {mutation.notation} applied successfully.")
        print(f"   Context: ...{sequence[max(0, idx-5):idx]}[{actual_aa}→{mutation.mutant_aa}]{sequence[idx+1:idx+6]}...")

        return mutated_sequence, mutation
