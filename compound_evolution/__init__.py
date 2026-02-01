"""
Compound Evolution Analyzer

A tool for analyzing chemical compound evolution in patents to identify
lead compounds and their optimization pathways.
"""

from .analyzer import CompoundEvolutionAnalyzer
from .fingerprints import FingerprintCalculator
from .clustering import CompoundClusterer
from .pathway import PathwayAnalyzer
from .visualization import EvolutionVisualizer

__version__ = "1.0.0"
__all__ = [
    "CompoundEvolutionAnalyzer",
    "FingerprintCalculator",
    "CompoundClusterer",
    "PathwayAnalyzer",
    "EvolutionVisualizer",
]
