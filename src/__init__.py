"""
Mutated Protein Modeling Pipeline for Binder Design

A complete pipeline for:
1. Fetching Wild Type structures from AlphaFold Database
2. Applying point mutations to sequences
3. Generating conformations using BioEmu
4. Performing structural analysis (Core vs Global RMSD)
5. Selecting the best candidate for RFdiffusion
"""

from .models import MutationInfo, CoreResidueInfo, RMSDResult
from .handlers import ProteinDataHandler
from .runners import BioEmuRunner
from .analysis import StructureAnalyzer
from .visualization import StructureVisualizer
from .exporters import MutantPDBExporter
from .pipeline import MutationPipeline

__version__ = "1.0.0"
__all__ = [
    "MutationInfo",
    "CoreResidueInfo", 
    "RMSDResult",
    "ProteinDataHandler",
    "BioEmuRunner",
    "StructureAnalyzer",
    "StructureVisualizer",
    "MutantPDBExporter",
    "MutationPipeline",
]
