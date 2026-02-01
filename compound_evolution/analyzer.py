"""
Main analyzer module for compound evolution analysis.

This module provides the high-level interface for analyzing chemical
compound evolution in patents. It orchestrates the parsing, fingerprinting,
clustering, pathway analysis, and visualization components.
"""

from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass, field
from pathlib import Path
import json
import numpy as np

from .fingerprints import Compound, FingerprintCalculator
from .clustering import Cluster, CompoundClusterer
from .pathway import EvolutionPathway, PathwayAnalyzer, LeadCompound, LeadType
from .visualization import EvolutionVisualizer


@dataclass
class AnalysisResult:
    """
    Complete results from compound evolution analysis.

    Attributes:
        compounds: All parsed compounds
        similarity_matrix: Pairwise Tanimoto similarities
        clusters: Identified compound clusters
        pathway: Evolution pathway with leads
        summary: Dictionary of summary statistics
    """
    compounds: List[Compound]
    similarity_matrix: np.ndarray
    clusters: List[Cluster]
    pathway: EvolutionPathway
    summary: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert results to a serializable dictionary."""
        return {
            'num_compounds': len(self.compounds),
            'num_clusters': len(self.clusters),
            'num_leads': len(self.pathway.all_leads),
            'initial_leads': [
                {
                    'name': l.compound.name,
                    'smiles': l.compound.smiles,
                    'generation': l.generation,
                    'analog_count': l.analog_count
                }
                for l in self.pathway.initial_leads
            ],
            'intermediate_leads': [
                {
                    'name': l.compound.name,
                    'smiles': l.compound.smiles,
                    'generation': l.generation,
                    'analog_count': l.analog_count
                }
                for l in self.pathway.intermediate_leads
            ],
            'final_lead': {
                'name': self.pathway.final_lead.compound.name,
                'smiles': self.pathway.final_lead.compound.smiles,
                'generation': self.pathway.final_lead.generation,
                'analog_count': self.pathway.final_lead.analog_count
            } if self.pathway.final_lead else None,
            'clusters': [
                {
                    'id': c.id,
                    'size': len(c.compounds),
                    'centroid': c.centroid.name if c.centroid else None,
                    'avg_similarity': float(c.avg_internal_similarity),
                    'density': float(c.density)
                }
                for c in self.clusters
            ],
            'summary': self.summary
        }

    def save_json(self, path: str) -> None:
        """Save results to a JSON file."""
        with open(path, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)


class CompoundEvolutionAnalyzer:
    """
    High-level analyzer for compound evolution in patents.

    This class provides a simple interface to analyze chemical compound
    evolution, identify lead compounds, and visualize the optimization
    pathway.

    Example usage:
        >>> analyzer = CompoundEvolutionAnalyzer()
        >>> result = analyzer.analyze_sdf("patent_compounds.sdf")
        >>> analyzer.visualize(result, "evolution_pathway.png")

    Attributes:
        similarity_threshold: Tanimoto threshold for "close analog" (0.7)
        fingerprint_radius: Morgan fingerprint radius (2 = ECFP4)
        fingerprint_bits: Number of fingerprint bits (2048)
        clustering_method: Clustering algorithm to use
    """

    CLUSTERING_METHODS = ['hierarchical', 'network', 'butina']

    def __init__(self,
                 similarity_threshold: float = 0.7,
                 fingerprint_radius: int = 2,
                 fingerprint_bits: int = 2048,
                 clustering_method: str = 'butina'):
        """
        Initialize the analyzer.

        Args:
            similarity_threshold: Minimum Tanimoto similarity for compounds
                                to be considered structurally related.
                                Higher values = stricter grouping.
            fingerprint_radius: Radius for Morgan fingerprints.
                              2 = ECFP4, 3 = ECFP6.
            fingerprint_bits: Length of fingerprint bit vectors.
            clustering_method: Algorithm for clustering compounds.
                             Options: 'hierarchical', 'network', 'butina'
        """
        if clustering_method not in self.CLUSTERING_METHODS:
            raise ValueError(
                f"Unknown clustering method: {clustering_method}. "
                f"Options: {self.CLUSTERING_METHODS}"
            )

        self.similarity_threshold = similarity_threshold
        self.fingerprint_radius = fingerprint_radius
        self.fingerprint_bits = fingerprint_bits
        self.clustering_method = clustering_method

        # Initialize components
        self.fp_calculator = FingerprintCalculator(
            radius=fingerprint_radius,
            n_bits=fingerprint_bits
        )
        self.clusterer = CompoundClusterer(
            similarity_threshold=similarity_threshold
        )
        self.pathway_analyzer = PathwayAnalyzer(
            similarity_threshold=similarity_threshold
        )
        self.visualizer = EvolutionVisualizer()

    def analyze_sdf(self, sdf_path: str) -> AnalysisResult:
        """
        Analyze compounds from an SDF file.

        This is the main entry point for analyzing a patent's compounds.
        It performs the complete analysis pipeline:
        1. Parse SDF and calculate fingerprints
        2. Compute similarity matrix
        3. Cluster compounds
        4. Identify lead compounds and pathway

        Args:
            sdf_path: Path to the SDF file containing compounds

        Returns:
            AnalysisResult with all analysis data

        Raises:
            FileNotFoundError: If SDF file doesn't exist
            ValueError: If no valid compounds found
        """
        # Step 1: Parse compounds
        print(f"Parsing compounds from {sdf_path}...")
        compounds = self.fp_calculator.parse_sdf(sdf_path)
        print(f"  Found {len(compounds)} valid compounds")

        return self._analyze_compounds(compounds)

    def analyze_smiles(self,
                      smiles_list: List[str],
                      names: Optional[List[str]] = None) -> AnalysisResult:
        """
        Analyze compounds from SMILES strings.

        Args:
            smiles_list: List of SMILES strings
            names: Optional list of compound names

        Returns:
            AnalysisResult with all analysis data
        """
        print(f"Parsing {len(smiles_list)} SMILES strings...")
        compounds = self.fp_calculator.parse_smiles_list(smiles_list, names)
        print(f"  Parsed {len(compounds)} valid compounds")

        return self._analyze_compounds(compounds)

    def _analyze_compounds(self, compounds: List[Compound]) -> AnalysisResult:
        """
        Perform analysis on parsed compounds.

        Args:
            compounds: List of Compound objects

        Returns:
            AnalysisResult with complete analysis
        """
        # Step 2: Calculate similarity matrix
        print("Calculating similarity matrix...")
        similarity_matrix = self.fp_calculator.calculate_similarity_matrix(
            compounds
        )

        # Step 3: Cluster compounds
        print(f"Clustering compounds using {self.clustering_method} method...")
        if self.clustering_method == 'hierarchical':
            clusters = self.clusterer.cluster_hierarchical(
                compounds, similarity_matrix
            )
        elif self.clustering_method == 'network':
            clusters = self.clusterer.cluster_network(
                compounds, similarity_matrix
            )
        else:  # butina
            clusters = self.clusterer.cluster_butina(
                compounds, similarity_matrix
            )
        print(f"  Identified {len(clusters)} clusters")

        # Step 4: Analyze pathway
        print("Analyzing evolution pathway...")
        pathway = self.pathway_analyzer.analyze(
            compounds, clusters, similarity_matrix
        )

        # Compile summary statistics
        summary = self._compute_summary(
            compounds, clusters, pathway, similarity_matrix
        )

        return AnalysisResult(
            compounds=compounds,
            similarity_matrix=similarity_matrix,
            clusters=clusters,
            pathway=pathway,
            summary=summary
        )

    def _compute_summary(self,
                        compounds: List[Compound],
                        clusters: List[Cluster],
                        pathway: EvolutionPathway,
                        similarity_matrix: np.ndarray) -> Dict[str, Any]:
        """
        Compute summary statistics for the analysis.

        Args:
            compounds: All compounds
            clusters: Compound clusters
            pathway: Evolution pathway
            similarity_matrix: Pairwise similarities

        Returns:
            Dictionary of summary statistics
        """
        # Overall similarity statistics
        triu_indices = np.triu_indices(len(compounds), k=1)
        pairwise_sims = similarity_matrix[triu_indices]

        summary = {
            'total_compounds': len(compounds),
            'total_clusters': len(clusters),
            'total_leads': len(pathway.all_leads),
            'avg_pairwise_similarity': float(np.mean(pairwise_sims)),
            'max_pairwise_similarity': float(np.max(pairwise_sims)),
            'min_pairwise_similarity': float(np.min(pairwise_sims)),
            'similarity_std': float(np.std(pairwise_sims)),
        }

        # Cluster statistics
        cluster_sizes = [len(c.compounds) for c in clusters]
        summary['avg_cluster_size'] = float(np.mean(cluster_sizes))
        summary['max_cluster_size'] = int(np.max(cluster_sizes))
        summary['min_cluster_size'] = int(np.min(cluster_sizes))

        # Pathway statistics
        if pathway.all_leads:
            summary['num_generations'] = max(
                l.generation for l in pathway.all_leads
            ) + 1
            summary['initial_lead_smiles'] = (
                pathway.initial_leads[0].compound.smiles
                if pathway.initial_leads else None
            )
            summary['final_lead_smiles'] = (
                pathway.final_lead.compound.smiles
                if pathway.final_lead else None
            )

        # Molecular property statistics
        mol_weights = [c.mol_weight for c in compounds]
        summary['avg_mol_weight'] = float(np.mean(mol_weights))
        summary['mol_weight_range'] = [
            float(np.min(mol_weights)),
            float(np.max(mol_weights))
        ]

        return summary

    def visualize(self,
                 result: AnalysisResult,
                 output_path: str,
                 style: str = 'network',
                 **kwargs) -> None:
        """
        Generate visualization of the analysis results.

        Args:
            result: Analysis results to visualize
            output_path: Path to save the visualization
            style: Visualization style ('network' or 'linear')
            **kwargs: Additional arguments for the visualizer

        Supported kwargs:
            show_structures: bool - Show chemical structures (default True)
            show_all_compounds: bool - Show all compounds (default True)
            title: str - Figure title
            figsize: tuple - Figure size (width, height)
            dpi: int - Output resolution
        """
        print(f"Generating {style} visualization...")

        if style == 'network':
            self.visualizer.visualize_pathway(
                pathway=result.pathway,
                compounds=result.compounds,
                similarity_matrix=result.similarity_matrix,
                clusters=result.clusters,
                output_path=output_path,
                **kwargs
            )
        elif style == 'linear':
            self.visualizer.visualize_pathway_linear(
                pathway=result.pathway,
                output_path=output_path,
                **kwargs
            )
        else:
            raise ValueError(f"Unknown visualization style: {style}")

        print(f"  Saved to {output_path}")

    def generate_report(self,
                       result: AnalysisResult,
                       output_dir: str,
                       base_name: str = "evolution_analysis") -> Dict[str, str]:
        """
        Generate comprehensive analysis report.

        Creates multiple output files:
        - Network visualization (PNG)
        - Linear pathway visualization (PNG)
        - HTML summary report
        - JSON data export

        Args:
            result: Analysis results
            output_dir: Directory for output files
            base_name: Base name for output files

        Returns:
            Dictionary mapping output type to file path
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        outputs = {}

        # Network visualization
        network_path = output_path / f"{base_name}_network.png"
        self.visualize(result, str(network_path), style='network')
        outputs['network_visualization'] = str(network_path)

        # Linear visualization
        linear_path = output_path / f"{base_name}_linear.png"
        self.visualize(result, str(linear_path), style='linear')
        outputs['linear_visualization'] = str(linear_path)

        # HTML report
        html_path = output_path / f"{base_name}_report.html"
        self.visualizer.create_summary_report(
            result.pathway,
            result.compounds,
            result.clusters,
            str(html_path)
        )
        outputs['html_report'] = str(html_path)

        # JSON export
        json_path = output_path / f"{base_name}_data.json"
        result.save_json(str(json_path))
        outputs['json_data'] = str(json_path)

        print(f"\nGenerated report files in {output_dir}:")
        for output_type, path in outputs.items():
            print(f"  - {output_type}: {Path(path).name}")

        return outputs

    def print_summary(self, result: AnalysisResult) -> None:
        """
        Print a text summary of the analysis results.

        Args:
            result: Analysis results
        """
        print("\n" + "=" * 60)
        print("COMPOUND EVOLUTION ANALYSIS SUMMARY")
        print("=" * 60)

        print(f"\nDataset Overview:")
        print(f"  Total compounds: {len(result.compounds)}")
        print(f"  Total clusters: {len(result.clusters)}")
        print(f"  Identified leads: {len(result.pathway.all_leads)}")

        print(f"\nSimilarity Statistics:")
        print(f"  Average pairwise similarity: "
              f"{result.summary['avg_pairwise_similarity']:.3f}")
        print(f"  Similarity range: "
              f"{result.summary['min_pairwise_similarity']:.3f} - "
              f"{result.summary['max_pairwise_similarity']:.3f}")

        print(f"\nCluster Statistics:")
        print(f"  Average cluster size: {result.summary['avg_cluster_size']:.1f}")
        print(f"  Cluster size range: "
              f"{result.summary['min_cluster_size']} - "
              f"{result.summary['max_cluster_size']}")

        print(f"\nEvolution Pathway:")
        print(f"  Generations identified: {result.summary.get('num_generations', 0)}")

        if result.pathway.initial_leads:
            print(f"\n  Initial Lead(s):")
            for lead in result.pathway.initial_leads:
                print(f"    - {lead.compound.name}")
                print(f"      SMILES: {lead.compound.smiles[:60]}...")
                print(f"      Analogs: {lead.analog_count}")

        if result.pathway.intermediate_leads:
            print(f"\n  Intermediate Lead(s):")
            for lead in result.pathway.intermediate_leads:
                print(f"    - {lead.compound.name} (Gen {lead.generation})")

        if result.pathway.final_lead:
            print(f"\n  Final Lead:")
            lead = result.pathway.final_lead
            print(f"    - {lead.compound.name}")
            print(f"      SMILES: {lead.compound.smiles[:60]}...")
            print(f"      Analogs: {lead.analog_count}")
            print(f"      Cluster density: {lead.cluster.density:.3f}")

        print("\n" + "=" * 60)
