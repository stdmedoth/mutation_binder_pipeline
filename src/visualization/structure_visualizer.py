"""
Protein structure visualizer using py3Dmol and matplotlib.
"""

from pathlib import Path
from typing import List, Optional, Tuple, Any

import numpy as np
import matplotlib.pyplot as plt
import py3Dmol

from ..models import RMSDResult


class StructureVisualizer:
    """
    Protein structure visualizer using py3Dmol.
    """

    @staticmethod
    def plot_pca_projection(
        pca_components: np.ndarray,
        sample_ids: List[str],
        pca_model: Any,
        wt_projected: Optional[np.ndarray] = None,
        figsize: Tuple[int, int] = (10, 8)
    ) -> plt.Figure:
        """
        Generate PCA scatter plot (Conformational Space).

        Args:
            pca_components: Matrix (N_samples, 2) with PC1/PC2 coordinates.
            sample_ids: List of sample IDs.
            pca_model: Trained PCA object (for explained variance).
            wt_projected: WT coordinates (1, 2) projected (optional).
            figsize: Figure size.

        Returns:
            Matplotlib Figure object.
        """
        fig, ax = plt.subplots(figsize=figsize)

        # 1. Plot Mutants (BioEmu samples)
        scatter = ax.scatter(
            pca_components[:, 0],
            pca_components[:, 1],
            c='blue',
            alpha=0.6,
            edgecolors='k',
            linewidth=0.5,
            s=80,
            label='BioEmu Conformations'
        )

        # 2. Plot WT (Reference) if provided
        if wt_projected is not None:
            ax.scatter(
                wt_projected[0, 0],
                wt_projected[0, 1],
                c='red',
                s=300,
                marker='*',
                edgecolors='black',
                linewidth=1.5,
                label='Wild Type (Ref)'
            )

            # Add WT annotation
            ax.annotate(
                "WT",
                (wt_projected[0, 0], wt_projected[0, 1]),
                xytext=(10, 10), textcoords='offset points',
                fontsize=12, fontweight='bold', color='red'
            )

        # 3. Labels and styling
        var_exp = pca_model.explained_variance_ratio_
        ax.set_xlabel(f'Principal Component 1 ({var_exp[0]*100:.1f}%)', fontsize=12)
        ax.set_ylabel(f'Principal Component 2 ({var_exp[1]*100:.1f}%)', fontsize=12)
        ax.set_title(
            f'Conformational Space Analysis (PCA)\nTotal Variance Explained: {sum(var_exp)*100:.1f}%',
            fontsize=14, fontweight='bold'
        )

        # Grid and legend
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(fontsize=11)

        # Add IDs to extreme points (outliers)
        std_x, std_y = np.std(pca_components[:, 0]), np.std(pca_components[:, 1])
        mean_x, mean_y = np.mean(pca_components[:, 0]), np.mean(pca_components[:, 1])

        for i, txt in enumerate(sample_ids):
            x, y = pca_components[i, 0], pca_components[i, 1]
            # Annotate if far from center
            if abs(x - mean_x) > 1.5 * std_x or abs(y - mean_y) > 1.5 * std_y:
                ax.annotate(txt, (x, y), fontsize=9, alpha=0.8, xytext=(5, 5), textcoords='offset points')

        plt.tight_layout()
        return fig

    @staticmethod
    def plot_rmsd_comparison(
        results: List[RMSDResult],
        best_id: Optional[str] = None,
        figsize: Tuple[int, int] = (14, 5)
    ) -> plt.Figure:
        """
        Generate bar charts comparing RMSD of conformations.

        Args:
            results: List of RMSDResults.
            best_id: ID of best conformation to highlight.
            figsize: Figure size.

        Returns:
            Matplotlib figure.
        """
        fig, axes = plt.subplots(1, 2, figsize=figsize)

        sample_ids = [r.sample_id for r in results]
        rmsd_core = [r.rmsd_core for r in results]
        rmsd_global = [r.rmsd_global for r in results]

        # Colors: highlight best candidate
        colors_core = ['#2ecc71' if r.sample_id == best_id else '#3498db' for r in results]
        colors_global = ['#2ecc71' if r.sample_id == best_id else '#e74c3c' for r in results]

        # Core RMSD plot
        ax1 = axes[0]
        ax1.bar(sample_ids, rmsd_core, color=colors_core, edgecolor='black', linewidth=0.5)
        ax1.set_xlabel('Conformation', fontsize=12)
        ax1.set_ylabel('Core RMSD (Å)', fontsize=12)
        ax1.set_title('Structured Core RMSD\n(pLDDT ≥ threshold)', fontsize=14, fontweight='bold')
        ax1.tick_params(axis='x', rotation=45)
        ax1.axhline(y=np.mean(rmsd_core), color='gray', linestyle='--', alpha=0.7, 
                    label=f'Mean: {np.mean(rmsd_core):.2f}Å')
        ax1.legend()

        # Global RMSD plot
        ax2 = axes[1]
        ax2.bar(sample_ids, rmsd_global, color=colors_global, edgecolor='black', linewidth=0.5)
        ax2.set_xlabel('Conformation', fontsize=12)
        ax2.set_ylabel('Global RMSD (Å)', fontsize=12)
        ax2.set_title('Global RMSD\n(after Core alignment)', fontsize=14, fontweight='bold')
        ax2.tick_params(axis='x', rotation=45)
        ax2.axhline(y=np.mean(rmsd_global), color='gray', linestyle='--', alpha=0.7, 
                    label=f'Mean: {np.mean(rmsd_global):.2f}Å')
        ax2.legend()

        plt.tight_layout()
        return fig

    @staticmethod
    def visualize_structures_3d(
        wt_pdb_path: Path,
        mutant_pdb_path: Path,
        mutation_position: int,
        width: int = 800,
        height: int = 600
    ):
        """
        Interactive 3D visualization with py3Dmol.

        - WT: gray, transparent
        - Mutant: green, solid
        - Mutated residue: red, stick

        Args:
            wt_pdb_path: Path to WT PDB.
            mutant_pdb_path: Path to mutant PDB.
            mutation_position: Mutation position (1-indexed).
            width: Viewer width.
            height: Viewer height.

        Returns:
            py3Dmol view object.
        """
        # Read PDB files
        with open(wt_pdb_path, 'r') as f:
            wt_pdb = f.read()
        with open(mutant_pdb_path, 'r') as f:
            mutant_pdb = f.read()

        # Create viewer
        view = py3Dmol.view(width=width, height=height)

        # Add WT (gray, transparent)
        view.addModel(wt_pdb, 'pdb')
        view.setStyle(
            {'model': 0},
            {'cartoon': {'color': 'gray', 'opacity': 0.5}}
        )

        # Add Mutant (green, solid)
        view.addModel(mutant_pdb, 'pdb')
        view.setStyle(
            {'model': 1},
            {'cartoon': {'color': 'green', 'opacity': 1.0}}
        )

        # Highlight mutated residue (red, stick)
        view.addStyle(
            {'model': 1, 'resi': mutation_position},
            {'stick': {'color': 'red', 'radius': 0.3}}
        )

        # Add label at mutated residue
        view.addLabel(
            f"Mutation (pos {mutation_position})",
            {'fontSize': 12, 'fontColor': 'white', 'backgroundColor': 'red'},
            {'model': 1, 'resi': mutation_position}
        )

        # Visualization settings
        view.zoomTo()
        view.setBackgroundColor('white')

        return view

    @staticmethod
    def plot_flexibility_analysis(
        results: List[RMSDResult],
        figsize: Tuple[int, int] = (10, 6)
    ) -> plt.Figure:
        """
        Flexibility analysis: scatter plot Core vs Global RMSD.

        Args:
            results: List of RMSDResults.
            figsize: Figure size.

        Returns:
            Matplotlib figure.
        """
        fig, ax = plt.subplots(figsize=figsize)

        rmsd_core = [r.rmsd_core for r in results]
        rmsd_global = [r.rmsd_global for r in results]
        sample_ids = [r.sample_id for r in results]

        # Scatter plot
        scatter = ax.scatter(
            rmsd_core, rmsd_global,
            c=rmsd_core, cmap='RdYlGn_r',
            s=100, edgecolor='black', linewidth=1
        )

        # Diagonal line (equal RMSD)
        max_val = max(max(rmsd_core), max(rmsd_global)) * 1.1
        ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='Core = Global')

        # Labels
        for i, sample_id in enumerate(sample_ids):
            ax.annotate(
                sample_id,
                (rmsd_core[i], rmsd_global[i]),
                textcoords="offset points",
                xytext=(5, 5),
                fontsize=8
            )

        ax.set_xlabel('Core RMSD (Å)', fontsize=12)
        ax.set_ylabel('Global RMSD (Å)', fontsize=12)
        ax.set_title('Flexibility Analysis: Core vs Global RMSD', fontsize=14, fontweight='bold')

        # Colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Core RMSD (Å)', fontsize=10)

        ax.legend()
        ax.set_xlim(0, max_val)
        ax.set_ylim(0, max_val)

        plt.tight_layout()
        return fig
