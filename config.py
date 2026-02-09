"""
Configuration settings for the Mutation Pipeline.
"""

from pathlib import Path

# ╔════════════════════════════════════════════════════════════════════╗
# ║                    DIRECTORY CONFIGURATION                         ║
# ╚════════════════════════════════════════════════════════════════════╝

# Base directory for all pipeline outputs
# Modify this to your preferred location
BASE_DIR = Path("./mutation_pipeline_results")

# Subdirectories (created automatically)
STRUCTURES_DIR = BASE_DIR / "structures"
BIOEMU_DIR = BASE_DIR / "bioemu_results"
OUTPUT_DIR = BASE_DIR / "output"


# ╔════════════════════════════════════════════════════════════════════╗
# ║                    PIPELINE PARAMETERS                             ║
# ╚════════════════════════════════════════════════════════════════════╝

# Number of conformations to generate with BioEmu
NUM_CONFORMATIONS = 10

# pLDDT threshold for Core detection
# 70: Includes moderately structured regions
# 80: Only well-structured regions
PLDDT_THRESHOLD = 70.0


# ╔════════════════════════════════════════════════════════════════════╗
# ║                    DEFAULT INPUT PARAMETERS                        ║
# ╚════════════════════════════════════════════════════════════════════╝

# UniProt ID of the protein of interest
# Examples: P04637 (p53), P38398 (BRCA1), P00533 (EGFR)
DEFAULT_UNIPROT_ID = "P50995"

# Point mutation in format: [original AA][position][mutant AA]
# Examples: R175H, G245S, R248Q (common p53 mutations)
DEFAULT_MUTATION = "P36R"


# ╔════════════════════════════════════════════════════════════════════╗
# ║                    VISUALIZATION SETTINGS                          ║
# ╚════════════════════════════════════════════════════════════════════╝

# Figure sizes for matplotlib plots
FIGURE_SIZE_RMSD = (14, 5)
FIGURE_SIZE_PCA = (10, 8)
FIGURE_SIZE_FLEXIBILITY = (10, 6)

# 3D viewer dimensions
VIEWER_WIDTH = 800
VIEWER_HEIGHT = 600

# Plot DPI for saved figures
PLOT_DPI = 150


# ╔════════════════════════════════════════════════════════════════════╗
# ║                    GOOGLE COLAB CONFIGURATION                      ║
# ╚════════════════════════════════════════════════════════════════════╝

# Set to True when running in Google Colab
USE_GOOGLE_DRIVE = False

# Google Drive paths (only used if USE_GOOGLE_DRIVE is True)
DRIVE_FOLDER = "/content/drive"
MY_DRIVE_PATH = f"{DRIVE_FOLDER}/MyDrive"
COLAB_FOLDER_NAME = "Projeto_Mutacoes"
COLAB_ROOT = Path(f"{MY_DRIVE_PATH}/{COLAB_FOLDER_NAME}")


def get_base_dir() -> Path:
    """
    Get the appropriate base directory based on environment.
    
    Returns:
        Path to base directory
    """
    if USE_GOOGLE_DRIVE:
        return COLAB_ROOT / "mutation_pipeline_results"
    return BASE_DIR


def setup_directories() -> None:
    """Create all necessary directories."""
    base = get_base_dir()
    for subdir in ["structures", "bioemu_results", "output"]:
        (base / subdir).mkdir(parents=True, exist_ok=True)
    print(f"📁 Directories created at: {base.absolute()}")
