"""
Tests for the Compound Evolution Analyzer.
"""

import pytest
import numpy as np
import tempfile
import os

from compound_evolution import (
    CompoundEvolutionAnalyzer,
    FingerprintCalculator,
    CompoundClusterer,
    PathwayAnalyzer,
)
from compound_evolution.fingerprints import Compound
from compound_evolution.pathway import LeadType


# Test data - small SAR series
TEST_SMILES = [
    ("c1ccc(Nc2ncnc3ccccc23)cc1", "Parent"),
    ("Fc1ccc(Nc2ncnc3ccccc23)cc1", "4-Fluoro"),
    ("Fc1cccc(Nc2ncnc3ccccc23)c1", "3-Fluoro"),
    ("Fc1ccc(Nc2ncnc3ccccc23)cc1F", "3,4-DiFluoro"),
    ("Clc1ccc(Nc2ncnc3ccccc23)cc1", "4-Chloro"),
    ("COc1ccc(Nc2ncnc3ccccc23)cc1", "4-Methoxy"),
    ("Cc1ccc(Nc2ncnc3ccccc23)cc1", "4-Methyl"),
    ("c1ccc(Nc2ncnc3cc(Cl)ccc23)cc1", "Core-6-Cl"),
]


class TestFingerprintCalculator:
    """Tests for FingerprintCalculator class."""

    def test_init_default_params(self):
        """Test initialization with default parameters."""
        calc = FingerprintCalculator()
        assert calc.radius == 2
        assert calc.n_bits == 2048

    def test_init_custom_params(self):
        """Test initialization with custom parameters."""
        calc = FingerprintCalculator(radius=3, n_bits=1024)
        assert calc.radius == 3
        assert calc.n_bits == 1024

    def test_parse_smiles_list(self):
        """Test parsing SMILES strings."""
        calc = FingerprintCalculator()
        smiles = [s for s, _ in TEST_SMILES]
        names = [n for _, n in TEST_SMILES]

        compounds = calc.parse_smiles_list(smiles, names)

        assert len(compounds) == len(TEST_SMILES)
        assert all(isinstance(c, Compound) for c in compounds)
        assert compounds[0].name == "Parent"

    def test_parse_smiles_invalid(self):
        """Test handling of invalid SMILES."""
        calc = FingerprintCalculator()
        smiles = ["c1ccccc1", "invalid_smiles", "CCO"]

        compounds = calc.parse_smiles_list(smiles)

        assert len(compounds) == 2  # Only valid ones

    def test_calculate_similarity(self):
        """Test Tanimoto similarity calculation."""
        calc = FingerprintCalculator()
        smiles = [s for s, _ in TEST_SMILES[:2]]
        compounds = calc.parse_smiles_list(smiles)

        sim = calc.calculate_similarity(compounds[0], compounds[1])

        assert 0.0 <= sim <= 1.0
        assert sim > 0.7  # These should be quite similar

    def test_calculate_similarity_identical(self):
        """Test similarity of identical compounds is 1.0."""
        calc = FingerprintCalculator()
        compounds = calc.parse_smiles_list(["CCO"])

        sim = calc.calculate_similarity(compounds[0], compounds[0])

        assert sim == 1.0

    def test_calculate_similarity_matrix(self):
        """Test similarity matrix calculation."""
        calc = FingerprintCalculator()
        smiles = [s for s, _ in TEST_SMILES]
        compounds = calc.parse_smiles_list(smiles)

        matrix = calc.calculate_similarity_matrix(compounds)

        assert matrix.shape == (len(compounds), len(compounds))
        assert np.all(np.diag(matrix) == 1.0)  # Diagonal is 1
        assert np.allclose(matrix, matrix.T)  # Symmetric

    def test_get_similar_compounds(self):
        """Test finding similar compounds."""
        calc = FingerprintCalculator()
        smiles = [s for s, _ in TEST_SMILES]
        names = [n for _, n in TEST_SMILES]
        compounds = calc.parse_smiles_list(smiles, names)

        similar = calc.get_similar_compounds(
            compounds[0], compounds, threshold=0.7
        )

        assert len(similar) > 0
        assert all(sim >= 0.7 for _, sim in similar)


class TestCompoundClusterer:
    """Tests for CompoundClusterer class."""

    @pytest.fixture
    def compounds_and_matrix(self):
        """Create test compounds and similarity matrix."""
        calc = FingerprintCalculator()
        smiles = [s for s, _ in TEST_SMILES]
        names = [n for _, n in TEST_SMILES]
        compounds = calc.parse_smiles_list(smiles, names)
        matrix = calc.calculate_similarity_matrix(compounds)
        return compounds, matrix

    def test_cluster_hierarchical(self, compounds_and_matrix):
        """Test hierarchical clustering."""
        compounds, matrix = compounds_and_matrix
        clusterer = CompoundClusterer(similarity_threshold=0.7)

        clusters = clusterer.cluster_hierarchical(compounds, matrix)

        assert len(clusters) > 0
        total_compounds = sum(len(c.compounds) for c in clusters)
        assert total_compounds == len(compounds)

    def test_cluster_network(self, compounds_and_matrix):
        """Test network-based clustering."""
        compounds, matrix = compounds_and_matrix
        clusterer = CompoundClusterer(similarity_threshold=0.7)

        clusters = clusterer.cluster_network(compounds, matrix)

        assert len(clusters) > 0
        total_compounds = sum(len(c.compounds) for c in clusters)
        assert total_compounds == len(compounds)

    def test_cluster_butina(self, compounds_and_matrix):
        """Test Butina clustering."""
        compounds, matrix = compounds_and_matrix
        clusterer = CompoundClusterer(similarity_threshold=0.7)

        clusters = clusterer.cluster_butina(compounds, matrix)

        assert len(clusters) > 0
        total_compounds = sum(len(c.compounds) for c in clusters)
        assert total_compounds == len(compounds)

    def test_cluster_has_centroid(self, compounds_and_matrix):
        """Test that clusters have centroids assigned."""
        compounds, matrix = compounds_and_matrix
        clusterer = CompoundClusterer(similarity_threshold=0.7)

        clusters = clusterer.cluster_butina(compounds, matrix)

        for cluster in clusters:
            assert cluster.centroid is not None
            assert cluster.centroid in cluster.compounds


class TestPathwayAnalyzer:
    """Tests for PathwayAnalyzer class."""

    @pytest.fixture
    def analysis_data(self):
        """Create test data for pathway analysis."""
        calc = FingerprintCalculator()
        smiles = [s for s, _ in TEST_SMILES]
        names = [n for _, n in TEST_SMILES]
        compounds = calc.parse_smiles_list(smiles, names)
        matrix = calc.calculate_similarity_matrix(compounds)

        clusterer = CompoundClusterer(similarity_threshold=0.7)
        clusters = clusterer.cluster_butina(compounds, matrix)

        return compounds, clusters, matrix

    def test_analyze_pathway(self, analysis_data):
        """Test pathway analysis."""
        compounds, clusters, matrix = analysis_data
        analyzer = PathwayAnalyzer(similarity_threshold=0.7)

        pathway = analyzer.analyze(compounds, clusters, matrix)

        assert pathway is not None
        assert len(pathway.all_leads) > 0

    def test_pathway_has_initial_lead(self, analysis_data):
        """Test that pathway identifies initial lead."""
        compounds, clusters, matrix = analysis_data
        analyzer = PathwayAnalyzer(similarity_threshold=0.7)

        pathway = analyzer.analyze(compounds, clusters, matrix)

        assert len(pathway.initial_leads) > 0
        assert pathway.initial_leads[0].lead_type == LeadType.INITIAL

    def test_pathway_generations(self, analysis_data):
        """Test that generations are assigned."""
        compounds, clusters, matrix = analysis_data
        analyzer = PathwayAnalyzer(similarity_threshold=0.7)

        pathway = analyzer.analyze(compounds, clusters, matrix)

        generations = [lead.generation for lead in pathway.all_leads]
        assert 0 in generations  # At least initial generation


class TestCompoundEvolutionAnalyzer:
    """Tests for the main analyzer class."""

    def test_init_default(self):
        """Test initialization with defaults."""
        analyzer = CompoundEvolutionAnalyzer()

        assert analyzer.similarity_threshold == 0.7
        assert analyzer.clustering_method == 'butina'

    def test_init_custom(self):
        """Test initialization with custom parameters."""
        analyzer = CompoundEvolutionAnalyzer(
            similarity_threshold=0.8,
            clustering_method='network'
        )

        assert analyzer.similarity_threshold == 0.8
        assert analyzer.clustering_method == 'network'

    def test_init_invalid_clustering(self):
        """Test that invalid clustering method raises error."""
        with pytest.raises(ValueError):
            CompoundEvolutionAnalyzer(clustering_method='invalid')

    def test_analyze_smiles(self):
        """Test analysis from SMILES list."""
        analyzer = CompoundEvolutionAnalyzer()
        smiles = [s for s, _ in TEST_SMILES]
        names = [n for _, n in TEST_SMILES]

        result = analyzer.analyze_smiles(smiles, names)

        assert result is not None
        assert len(result.compounds) == len(TEST_SMILES)
        assert len(result.clusters) > 0
        assert result.pathway is not None

    def test_analyze_sdf(self):
        """Test analysis from SDF file."""
        from rdkit import Chem
        from rdkit.Chem import AllChem

        # Create temporary SDF file
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.sdf', delete=False
        ) as f:
            writer = Chem.SDWriter(f.name)
            for smiles, name in TEST_SMILES:
                mol = Chem.MolFromSmiles(smiles)
                AllChem.Compute2DCoords(mol)
                mol.SetProp("_Name", name)
                writer.write(mol)
            writer.close()
            temp_path = f.name

        try:
            analyzer = CompoundEvolutionAnalyzer()
            result = analyzer.analyze_sdf(temp_path)

            assert result is not None
            assert len(result.compounds) == len(TEST_SMILES)
        finally:
            os.unlink(temp_path)

    def test_result_to_dict(self):
        """Test converting result to dictionary."""
        analyzer = CompoundEvolutionAnalyzer()
        smiles = [s for s, _ in TEST_SMILES]
        names = [n for _, n in TEST_SMILES]

        result = analyzer.analyze_smiles(smiles, names)
        result_dict = result.to_dict()

        assert 'num_compounds' in result_dict
        assert 'num_clusters' in result_dict
        assert 'initial_leads' in result_dict
        assert 'final_lead' in result_dict

    def test_result_save_json(self):
        """Test saving result to JSON."""
        analyzer = CompoundEvolutionAnalyzer()
        smiles = [s for s, _ in TEST_SMILES]
        names = [n for _, n in TEST_SMILES]

        result = analyzer.analyze_smiles(smiles, names)

        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.json', delete=False
        ) as f:
            temp_path = f.name

        try:
            result.save_json(temp_path)
            assert os.path.exists(temp_path)
            assert os.path.getsize(temp_path) > 0
        finally:
            os.unlink(temp_path)


class TestVisualization:
    """Tests for visualization (basic smoke tests)."""

    @pytest.fixture
    def analysis_result(self):
        """Create analysis result for visualization tests."""
        analyzer = CompoundEvolutionAnalyzer()
        smiles = [s for s, _ in TEST_SMILES]
        names = [n for _, n in TEST_SMILES]
        return analyzer.analyze_smiles(smiles, names), analyzer

    def test_visualize_network(self, analysis_result):
        """Test network visualization generation."""
        result, analyzer = analysis_result

        with tempfile.NamedTemporaryFile(
            suffix='.png', delete=False
        ) as f:
            temp_path = f.name

        try:
            analyzer.visualize(result, temp_path, style='network')
            assert os.path.exists(temp_path)
            assert os.path.getsize(temp_path) > 0
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)

    def test_visualize_linear(self, analysis_result):
        """Test linear visualization generation."""
        result, analyzer = analysis_result

        with tempfile.NamedTemporaryFile(
            suffix='.png', delete=False
        ) as f:
            temp_path = f.name

        try:
            analyzer.visualize(result, temp_path, style='linear')
            assert os.path.exists(temp_path)
            assert os.path.getsize(temp_path) > 0
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)
