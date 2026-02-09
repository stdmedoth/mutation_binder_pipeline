"""
Complete pipeline for mutated protein modeling.
"""

from pathlib import Path
from typing import Dict, List, Optional

from ..handlers import ProteinDataHandler
from ..runners import BioEmuRunner
from ..analysis import StructureAnalyzer
from ..visualization import StructureVisualizer
from ..exporters import MutantPDBExporter
from ..models import MutationInfo, RMSDResult


class MutationPipeline:
    """
    Complete pipeline for mutated protein modeling.

    Orchestrates all stages:
    1. Data retrieval (sequence + WT structure)
    2. Mutation application
    3. Conformation generation (BioEmu)
    4. Structural analysis (RMSD)
    5. Selection and visualization
    6. Export for RFdiffusion
    """

    def __init__(
        self,
        base_dir: Path,
        plddt_threshold: float = 70.0,
        num_conformations: int = 10
    ):
        """
        Initialize the pipeline.

        Args:
            base_dir: Base directory for all outputs.
            plddt_threshold: pLDDT threshold for Core detection.
            num_conformations: Number of conformations to generate.
        """
        self.base_dir = Path(base_dir)
        self.plddt_threshold = plddt_threshold
        self.num_conformations = num_conformations

        # Create directories
        self.structures_dir = self.base_dir / "structures"
        self.bioemu_dir = self.base_dir / "bioemu_results"
        self.output_dir = self.base_dir / "output"
        
        for directory in [self.structures_dir, self.bioemu_dir, self.output_dir]:
            directory.mkdir(parents=True, exist_ok=True)

        # Initialize components
        self.data_handler = ProteinDataHandler(output_dir=self.structures_dir)
        self.bioemu_runner = BioEmuRunner(output_dir=self.bioemu_dir)
        self.analyzer = StructureAnalyzer(plddt_threshold=plddt_threshold)
        self.visualizer = StructureVisualizer()
        self.exporter = MutantPDBExporter()

        # Results
        self.wt_sequence: Optional[str] = None
        self.mutant_sequence: Optional[str] = None
        self.mutation_info: Optional[MutationInfo] = None
        self.wt_structure_path: Optional[Path] = None
        self.wt_conformations: List[Path] = []
        self.conformations: List[Path] = []
        self.rmsd_results: List[RMSDResult] = []
        self.best_result: Optional[RMSDResult] = None
        self.principal_components = None
        self.sample_ids = None
        self.pca_model = None

    def run(
        self,
        uniprot_id: str,
        mutation_str: str
    ) -> Dict:
        """
        Execute the complete pipeline.

        Args:
            uniprot_id: UniProt identifier of the protein.
            mutation_str: Mutation string (e.g., 'D40G').

        Returns:
            Dictionary with all results.
        """
        print("=" * 70)
        print("🧬 MUTATED PROTEIN MODELING PIPELINE")
        print("=" * 70)

        # Stage 1: Fetch data
        print("\n📥 STAGE 1: Fetching data...")
        print("-" * 40)

        self.wt_sequence = self.data_handler.fetch_sequence(uniprot_id)
        self.wt_structure_path = self.data_handler.fetch_alphafold_structure(uniprot_id)

        # Stage 2: Apply mutation
        print("\n🔧 STAGE 2: Applying mutation...")
        print("-" * 40)

        self.mutant_sequence, self.mutation_info = self.data_handler.apply_mutation(
            self.wt_sequence, mutation_str
        )

        # Stage 3: Generate conformations
        print("\n🔬 STAGE 3: Generating conformations with BioEmu...")
        print("-" * 40)

        self.wt_conformations = self.bioemu_runner.run_sampling(
            self.wt_sequence,
            num_samples=self.num_conformations,
            prefix=f"{uniprot_id}_WT"
        )

        self.conformations = self.bioemu_runner.run_sampling(
            self.mutant_sequence,
            num_samples=self.num_conformations,
            prefix=f"{uniprot_id}_{self.mutation_info.notation}"
        )

        # Stage 4: Structural analysis
        print("\n📊 STAGE 4: Structural analysis...")
        print("-" * 40)

        self.analyzer.load_wt_structure(self.wt_structure_path)
        self.analyzer.detect_core_residues()
        self.rmsd_results = self.analyzer.analyze_all_samples(
            self.wt_conformations, 
            self.conformations
        )
        self.best_result = self.analyzer.get_best_candidate(self.rmsd_results)

        self.principal_components, self.sample_ids, self.pca_model = \
            self.analyzer.calculate_pca(self.conformations)

        # Stage 5: Export result
        print("\n💾 STAGE 5: Exporting PDB for RFdiffusion...")
        print("-" * 40)

        output_pdb = self.output_dir / f"{uniprot_id}_{self.mutation_info.notation}_best_for_rfdiffusion.pdb"
        self.exporter.export_clean_pdb(self.best_result.sample_path, output_pdb)

        # Summary
        summary = self.exporter.generate_summary(
            uniprot_id=uniprot_id,
            mutation=self.mutation_info,
            best_result=self.best_result,
            all_results=self.rmsd_results,
            core_info=self.analyzer.core_info
        )
        print(summary)

        return {
            'uniprot_id': uniprot_id,
            'mutation': self.mutation_info,
            'wt_sequence': self.wt_sequence,
            'mutant_sequence': self.mutant_sequence,
            'wt_structure_path': self.wt_structure_path,
            'conformations': self.conformations,
            'rmsd_results': self.rmsd_results,
            'best_result': self.best_result,
            'output_pdb': output_pdb,
            'core_info': self.analyzer.core_info
        }
