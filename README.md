# 🧬 Mutation Binder Pipeline

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Google Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/)

Automated pipeline for mutated protein modeling and conformation selection for binder design using **RFdiffusion**.

## 📋 Table of Contents

- [Overview](#-overview)
- [Pipeline Architecture](#-pipeline-architecture)
- [Requirements](#-requirements)
- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Detailed Documentation](#-detailed-documentation)
- [Outputs](#-outputs)
- [Scientific Background](#-scientific-background)
- [Contributing](#-contributing)
- [License](#-license)

## 🎯 Overview

This pipeline was developed to automate the process of:

1. **Data retrieval** → Fetch sequences (UniProt) and Wild Type structures (AlphaFold DB)
2. **Mutation application** → Validation and application of point mutations
3. **Conformational sampling** → Generation of multiple conformations using BioEmu
4. **Intelligent structural analysis** → RMSD calculation with automatic Core detection
5. **Best candidate selection** → Identification of the ideal conformation for RFdiffusion

### Why this pipeline?

Mutated proteins often exhibit conformational changes that can affect their function and interactions. For binder design (protein ligands), it is crucial to:

- Obtain a representative structure of the mutated protein
- Ensure that the **central structural fold** is preserved
- Ignore natural fluctuations of disordered regions (IDRs)

## 🏗 Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                         INPUT                                        │
│              UniProt ID + Mutation (e.g., P04637, R175H)            │
└─────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────┐
│                    ProteinDataHandler                                │
│  ├── fetch_sequence()      → UniProt REST API                       │
│  ├── fetch_alphafold_structure() → AlphaFold DB v4                  │
│  └── apply_mutation()      → Validation + Application               │
└─────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────┐
│                       BioEmuRunner                                   │
│  └── run_sampling()        → N conformations (GPU required)         │
└─────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────┐
│                    StructureAnalyzer                                 │
│  ├── detect_core_residues() → Auto-Core via pLDDT                   │
│  ├── calculate_rmsd()       → Core + Global RMSD                    │
│  └── get_best_candidate()   → Selection by lowest Core RMSD        │
└─────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────┐
│                         OUTPUT                                       │
│  ├── Clean PDB for RFdiffusion                                      │
│  ├── Analysis plots (RMSD, flexibility)                             │
│  └── Interactive 3D visualization (py3Dmol)                         │
└─────────────────────────────────────────────────────────────────────┘
```

## 💻 Requirements

### Hardware
- **NVIDIA GPU** (required for BioEmu)
- Minimum 8GB VRAM recommended

### Software
- Python 3.10+
- CUDA 11.x or higher

### Python Libraries
```
biopython>=1.80
bioemu
py3dmol
matplotlib
seaborn
requests
numpy
```

## 🚀 Installation

### Option 1: Google Colab (Recommended)

1. Open the notebook in Google Colab
2. Set runtime to GPU: `Runtime > Change runtime type > GPU`
3. Run the installation cell

### Option 2: Local Installation

```bash
# Clone the repository
git clone https://github.com/your-username/mutation_binder_pipeline.git
cd mutation_binder_pipeline

# Create a virtual environment
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or: venv\Scripts\activate  # Windows

# Install dependencies
pip install biopython requests py3dmol matplotlib seaborn
pip install bioemu
```

## 📖 Quick Start

```python
from mutation_binder_pipeline import MutationPipeline

# Initialize pipeline
pipeline = MutationPipeline(
    plddt_threshold=70.0,    # Threshold for Core detection
    num_conformations=10      # Number of conformations to generate
)

# Run
results = pipeline.run(
    uniprot_id="P04637",      # p53 tumor suppressor
    mutation_str="R175H"       # Common cancer mutation
)

# Access results
print(f"Best conformation: {results['best_result'].sample_id}")
print(f"Core RMSD: {results['best_result'].rmsd_core:.3f} Å")
print(f"PDB for RFdiffusion: {results['output_pdb']}")
```

## 📚 Detailed Documentation

### Main Classes

#### `ProteinDataHandler`
Manages protein data retrieval and manipulation.

| Method | Description |
|--------|-------------|
| `fetch_sequence(uniprot_id)` | Fetches FASTA sequence from UniProt |
| `fetch_alphafold_structure(uniprot_id)` | Downloads PDB structure from AlphaFold DB |
| `apply_mutation(sequence, mutation_str)` | Validates and applies point mutation |

#### `BioEmuRunner`
Wrapper for the BioEmu generative model.

| Method | Description |
|--------|-------------|
| `initialize_model()` | Loads BioEmu model on GPU |
| `run_sampling(sequence, num_samples)` | Generates N structural conformations |

#### `StructureAnalyzer`
Structural analysis with intelligent Core detection.

| Method | Description |
|--------|-------------|
| `detect_core_residues()` | Identifies structured residues via pLDDT |
| `calculate_rmsd(sample_path)` | Calculates Core and Global RMSD |
| `get_best_candidate(results)` | Selects conformation with lowest Core RMSD |

### Mutation Format

The mutation string must follow the standard format:

```
[Original Amino Acid][Position][Mutant Amino Acid]
```

Examples:
- `R175H` → Arginine at position 175 → Histidine
- `D40G` → Aspartic Acid at position 40 → Glycine
- `G245S` → Glycine at position 245 → Serine

## 📊 Outputs

### Generated Files

```
mutation_pipeline_results/
├── structures/
│   └── AF-{UNIPROT_ID}-F1-model_v4.pdb     # WT structure (AlphaFold)
├── bioemu_results/
│   ├── {ID}_{MUTATION}_conf_000.pdb        # Conformation 1
│   ├── {ID}_{MUTATION}_conf_001.pdb        # Conformation 2
│   └── ...
└── output/
    ├── {ID}_{MUTATION}_best_for_rfdiffusion.pdb  # ⭐ Final PDB
    ├── {ID}_{MUTATION}_rmsd_comparison.png       # RMSD plot
    └── {ID}_{MUTATION}_flexibility_analysis.png  # Flexibility analysis
```

### 3D Visualization

The notebook generates an interactive visualization with py3Dmol:

| Element | Color | Style |
|---------|-------|-------|
| Wild Type | Gray | Cartoon, transparent |
| Mutant (best) | Green | Cartoon, solid |
| Mutated residue | Red | Stick |

## 🔬 Scientific Background

### Auto-Core Detection

AlphaFold stores the **pLDDT** (predicted Local Distance Difference Test) in the B-factor column of the PDB file. This score (0-100) indicates prediction confidence:

| pLDDT | Interpretation |
|-------|----------------|
| > 90 | Very high confidence |
| 70-90 | High confidence |
| 50-70 | Low confidence |
| < 50 | Very low (likely IDR) |

The pipeline uses pLDDT ≥ 70 (configurable) to define the **structured Core**.

### RMSD Calculation

RMSD (Root Mean Square Deviation) is calculated as:

$$\text{RMSD} = \sqrt{\frac{1}{N}\sum_{i=1}^{N}\|\vec{r}_i^{(A)} - \vec{r}_i^{(B)}\|^2}$$

Where:
- $\vec{r}_i^{(A)}$ and $\vec{r}_i^{(B)}$ are positions of atom $i$ in structures A and B
- Alignment is performed using only Cα atoms from the Core

### Why Separate Core and Global?

| Metric | What it measures | Importance for RFdiffusion |
|--------|------------------|----------------------------|
| **Core RMSD** | Deviation of structured fold | ⭐ Critical - should be minimal |
| **Global RMSD** | Total deviation including IDRs | Informative - expected to be higher |

## 🔗 Integration with RFdiffusion

The exported PDB is ready for use in RFdiffusion:

```bash
python scripts/run_inference.py \
    inference.input_pdb=output/P04637_R175H_best_for_rfdiffusion.pdb \
    inference.output_prefix=binder_design \
    'contigmap.contigs=[A1-100/0 B1-80]' \
    inference.num_designs=10
```

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the project
2. Create a branch for your feature (`git checkout -b feature/NewFeature`)
3. Commit your changes (`git commit -m 'Add: New feature'`)
4. Push to the branch (`git push origin feature/NewFeature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Contact

For questions or suggestions, open an [Issue](https://github.com/your-username/mutation_binder_pipeline/issues).

---

<p align="center">
  Developed with ❤️ for the Structural Bioinformatics community
</p>
