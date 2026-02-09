#!/usr/bin/env python3
"""
Mutated Protein Modeling Pipeline - Main Entry Point

This script runs the complete pipeline for:
1. Fetching Wild Type structures from AlphaFold Database
2. Applying point mutations to sequences
3. Generating conformations using BioEmu
4. Performing structural analysis (Core vs Global RMSD)
5. Selecting the best candidate for RFdiffusion

Usage:
    python main.py
    python main.py --uniprot P50995 --mutation P36R
    python main.py --uniprot P04637 --mutation R175H --conformations 20
"""

import argparse
import subprocess
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Suppress warnings
warnings.filterwarnings('ignore')

# Import configuration
from config import (
    DEFAULT_UNIPROT_ID,
    DEFAULT_MUTATION,
    NUM_CONFORMATIONS,
    PLDDT_THRESHOLD,
    PLOT_DPI,
    get_base_dir,
    setup_directories,
)


def check_gpu_availability() -> bool:
    """Check if GPU is available (required for BioEmu)."""
    try:
        result = subprocess.run(['nvidia-smi'], capture_output=True, text=True)
        if result.returncode == 0:
            print("✅ GPU detected!")
            print(result.stdout)
            return True
        else:
            print("❌ No GPU detected.")
            return False
    except FileNotFoundError:
        print("❌ nvidia-smi not found. Make sure you're in a GPU-enabled environment.")
        return False


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Mutated Protein Modeling Pipeline for Binder Design",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python main.py
    python main.py --uniprot P50995 --mutation P36R
    python main.py --uniprot P04637 --mutation R175H --conformations 20
    python main.py --uniprot P38398 --mutation C61G --plddt 80.0
        """
    )
    
    parser.add_argument(
        '--uniprot', '-u',
        type=str,
        default=DEFAULT_UNIPROT_ID,
        help=f'UniProt ID of the protein (default: {DEFAULT_UNIPROT_ID})'
    )
    
    parser.add_argument(
        '--mutation', '-m',
        type=str,
        default=DEFAULT_MUTATION,
        help=f'Mutation in format [AA][pos][AA] e.g., D40G (default: {DEFAULT_MUTATION})'
    )
    
    parser.add_argument(
        '--conformations', '-n',
        type=int,
        default=NUM_CONFORMATIONS,
        help=f'Number of conformations to generate (default: {NUM_CONFORMATIONS})'
    )
    
    parser.add_argument(
        '--plddt', '-p',
        type=float,
        default=PLDDT_THRESHOLD,
        help=f'pLDDT threshold for Core detection (default: {PLDDT_THRESHOLD})'
    )
    
    parser.add_argument(
        '--output-dir', '-o',
        type=str,
        default=None,
        help='Custom output directory (default: ./mutation_pipeline_results)'
    )
    
    parser.add_argument(
        '--skip-gpu-check',
        action='store_true',
        help='Skip GPU availability check'
    )
    
    parser.add_argument(
        '--no-plots',
        action='store_true',
        help='Skip generating visualization plots'
    )
    
    return parser.parse_args()


def run_visualizations(pipeline, output_dir: Path, uniprot_id: str, mutation: str):
    """Generate and save all visualization plots."""
    from src import StructureVisualizer
    
    print("\n📊 Generating visualizations...")
    print("-" * 40)
    
    # Set style
    sns.set_theme(style="whitegrid")
    plt.rcParams['figure.figsize'] = (12, 6)
    plt.rcParams['font.size'] = 12
    
    # 1. RMSD Comparison
    if pipeline.rmsd_results:
        fig = StructureVisualizer.plot_rmsd_comparison(
            pipeline.rmsd_results,
            best_id=pipeline.best_result.sample_id
        )
        fig.savefig(
            output_dir / f"{uniprot_id}_{mutation}_rmsd_comparison.png",
            dpi=PLOT_DPI, bbox_inches='tight'
        )
        plt.close(fig)
        print(f"   ✅ RMSD comparison plot saved")
    
    # 2. Flexibility Analysis
    if pipeline.rmsd_results:
        fig = StructureVisualizer.plot_flexibility_analysis(pipeline.rmsd_results)
        fig.savefig(
            output_dir / f"{uniprot_id}_{mutation}_flexibility_analysis.png",
            dpi=PLOT_DPI, bbox_inches='tight'
        )
        plt.close(fig)
        print(f"   ✅ Flexibility analysis plot saved")
    
    # 3. PCA Projection
    if pipeline.principal_components is not None:
        try:
            # Project WT
            wt_structure = pipeline.analyzer.wt_structure
            wt_all_cas = pipeline.analyzer._get_ca_atoms(wt_structure)
            wt_coords_flat = np.array([atom.get_coord() for atom in wt_all_cas]).flatten()
            wt_projected = pipeline.pca_model.transform([wt_coords_flat])
            
            fig = StructureVisualizer.plot_pca_projection(
                pca_components=pipeline.principal_components,
                sample_ids=pipeline.sample_ids,
                pca_model=pipeline.pca_model,
                wt_projected=wt_projected,
                figsize=(10, 8)
            )
            fig.savefig(
                output_dir / f"{uniprot_id}_{mutation}_pca_projection.png",
                dpi=PLOT_DPI, bbox_inches='tight'
            )
            plt.close(fig)
            print(f"   ✅ PCA projection plot saved")
        except Exception as e:
            print(f"   ⚠️ Could not generate PCA plot: {e}")
    
    print(f"\n📁 All plots saved to: {output_dir}")


def list_generated_files(base_dir: Path):
    """List all generated files."""
    print("\n📁 GENERATED FILES:")
    print("=" * 50)

    for subdir in ["structures", "bioemu_results", "output"]:
        directory = base_dir / subdir
        if directory.exists():
            files = list(directory.iterdir())
            if files:
                print(f"\n📂 {directory}:")
                for f in sorted(files):
                    if f.is_file():
                        size = f.stat().st_size / 1024  # KB
                        print(f"   └─ {f.name} ({size:.1f} KB)")
                    elif f.is_dir():
                        subfiles = list(f.iterdir())
                        print(f"   └─ {f.name}/ ({len(subfiles)} files)")


def main():
    """Main entry point."""
    args = parse_arguments()
    
    print("=" * 70)
    print("🧬 MUTATED PROTEIN MODELING PIPELINE")
    print("=" * 70)
    
    # Configuration summary
    print(f"\n📋 Configuration:")
    print(f"   UniProt ID: {args.uniprot}")
    print(f"   Mutation: {args.mutation}")
    print(f"   Conformations: {args.conformations}")
    print(f"   pLDDT Threshold: {args.plddt}")
    
    # GPU Check
    if not args.skip_gpu_check:
        print("\n🔍 Checking GPU availability...")
        gpu_available = check_gpu_availability()
        
        if not gpu_available:
            print("\n⚠️  WARNING: BioEmu requires GPU for efficient execution.")
            print("   Use --skip-gpu-check to bypass this check.")
            response = input("   Continue anyway? [y/N]: ")
            if response.lower() != 'y':
                print("   Exiting.")
                sys.exit(1)
    
    # Setup directories
    if args.output_dir:
        base_dir = Path(args.output_dir)
    else:
        base_dir = get_base_dir()
    
    for subdir in ["structures", "bioemu_results", "output"]:
        (base_dir / subdir).mkdir(parents=True, exist_ok=True)
    
    print(f"\n📁 Output directory: {base_dir.absolute()}")
    
    # Import pipeline (after GPU check to avoid import errors)
    from src import MutationPipeline
    
    # Run pipeline
    pipeline = MutationPipeline(
        base_dir=base_dir,
        plddt_threshold=args.plddt,
        num_conformations=args.conformations
    )
    
    results = pipeline.run(args.uniprot, args.mutation)
    
    # Generate visualizations
    if not args.no_plots:
        run_visualizations(
            pipeline,
            base_dir / "output",
            args.uniprot,
            args.mutation
        )
    
    # List files
    list_generated_files(base_dir)
    
    # Final message
    print("\n" + "=" * 70)
    print("✅ PIPELINE COMPLETED SUCCESSFULLY!")
    print("=" * 70)
    print(f"\n🎯 Best candidate for RFdiffusion:")
    print(f"   {results['output_pdb']}")
    print("\n📌 Next step: Use this PDB as input for RFdiffusion binder design.")
    
    return results


if __name__ == "__main__":
    main()
