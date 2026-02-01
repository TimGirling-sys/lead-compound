# Compound Evolution Analyzer

A Python tool for analyzing chemical compound evolution in patents to identify lead compounds and their optimization pathways.

## Overview

In drug discovery patents, researchers typically:
1. Start with an initial lead compound (a "hit" from screening)
2. Create structural analogs around it (building a "moat" of similar structures)
3. Test these analogs and identify improved compounds
4. Repeat this process iteratively until reaching a final optimized lead

This tool analyzes the chemical structures claimed in a patent to reconstruct this optimization journey, identifying:
- **Initial lead compound(s)** - the starting scaffold(s)
- **Intermediate leads** - key optimization milestones
- **Final lead compound** - the optimized result
- **Analog clusters** - groups of similar structures around each lead

## Features

- **SDF Parsing**: Read chemical structures from standard SDF files
- **Fingerprint Calculation**: Morgan/ECFP fingerprints for structural comparison
- **Similarity Analysis**: Tanimoto similarity calculations between all compounds
- **Multiple Clustering Algorithms**:
  - Butina (sphere exclusion) - best for identifying lead series
  - Hierarchical clustering - traditional approach
  - Network-based (Louvain community detection) - for complex datasets
- **Pathway Analysis**: Automatic identification of the optimization pathway
- **Publication-Quality Visualizations**:
  - Network view with clusters and pathway arrows
  - Linear pathway diagram
  - Chemical structure rendering
  - HTML summary reports

## Installation

### Prerequisites

- Python 3.8 or higher
- RDKit (for cheminformatics)

### Install via pip

```bash
# Clone the repository
git clone https://github.com/example/compound-evolution.git
cd compound-evolution

# Install the package
pip install -e .
```

### Install dependencies manually

```bash
# Using conda (recommended for RDKit)
conda create -n compound-evolution python=3.10
conda activate compound-evolution
conda install -c conda-forge rdkit

# Install other dependencies
pip install numpy scipy networkx matplotlib pillow
```

### Development installation

```bash
pip install -e ".[dev]"
```

## Quick Start

### Command Line Interface

```bash
# Analyze an SDF file and generate a report
compound-evolution analyze patent_compounds.sdf -o results/

# With custom parameters
compound-evolution analyze compounds.sdf -o output/ \
    --threshold 0.8 \
    --clustering network \
    --name my_analysis
```

### Python API

```python
from compound_evolution import CompoundEvolutionAnalyzer

# Initialize analyzer
analyzer = CompoundEvolutionAnalyzer(
    similarity_threshold=0.7,  # Tanimoto threshold for "close analogs"
    clustering_method='butina'  # Best for lead identification
)

# Analyze an SDF file
result = analyzer.analyze_sdf("patent_compounds.sdf")

# Print summary
analyzer.print_summary(result)

# Generate visualizations
analyzer.visualize(result, "pathway_network.png", style='network')
analyzer.visualize(result, "pathway_linear.png", style='linear')

# Generate full report
outputs = analyzer.generate_report(result, "output/", base_name="analysis")
```

### Analyze from SMILES

```python
from compound_evolution import CompoundEvolutionAnalyzer

smiles_list = [
    "c1ccc(Nc2ncnc3ccccc23)cc1",      # Parent
    "Fc1ccc(Nc2ncnc3ccccc23)cc1",     # 4-Fluoro analog
    "Clc1ccc(Nc2ncnc3ccccc23)cc1",    # 4-Chloro analog
    # ... more compounds
]

names = ["Parent", "4-Fluoro", "4-Chloro"]

analyzer = CompoundEvolutionAnalyzer()
result = analyzer.analyze_smiles(smiles_list, names)
```

## Example Output

### Network Visualization

The network visualization shows:
- **Nodes**: Each compound (larger nodes are lead compounds)
- **Colors**:
  - Green: Initial lead
  - Blue: Intermediate leads
  - Red: Final lead
  - Gray: Analogs
- **Edges**: Structural similarity connections
- **Cluster hulls**: Dashed outlines around related compounds
- **Purple arrows**: The optimization pathway

### Linear Pathway

The linear visualization shows the progression from initial to final lead with chemical structures and connecting arrows.

## Detailed Usage

### Similarity Threshold

The `similarity_threshold` parameter (default: 0.7) controls what constitutes a "close analog":
- **0.7**: Standard threshold - captures typical SAR analogs
- **0.8**: Stricter - only very close analogs grouped together
- **0.6**: Looser - may group more diverse structures

### Clustering Methods

1. **Butina** (default):
   - Best for identifying lead series
   - Cluster centers naturally identify lead compounds
   - Good for typical patent datasets

2. **Hierarchical**:
   - Traditional agglomerative clustering
   - Good when you want control over cluster hierarchy
   - Uses average linkage

3. **Network**:
   - Uses Louvain community detection
   - Good for large, complex datasets
   - Can identify non-spherical clusters

### Fingerprint Parameters

- **radius** (default: 2): Morgan fingerprint radius
  - 2 = ECFP4 (captures up to 4-bond environments)
  - 3 = ECFP6 (more detailed, may be better for similar structures)
- **n_bits** (default: 2048): Fingerprint length
  - Higher values reduce fingerprint collisions

## Output Files

When using `generate_report()`, the following files are created:

| File | Description |
|------|-------------|
| `*_network.png` | Network visualization with clusters and pathway |
| `*_linear.png` | Linear pathway diagram with structures |
| `*_report.html` | Interactive HTML summary report |
| `*_data.json` | Machine-readable analysis data |

## API Reference

### CompoundEvolutionAnalyzer

Main class for compound evolution analysis.

```python
analyzer = CompoundEvolutionAnalyzer(
    similarity_threshold=0.7,    # Tanimoto cutoff for analogs
    fingerprint_radius=2,        # Morgan fingerprint radius
    fingerprint_bits=2048,       # Fingerprint length
    clustering_method='butina'   # 'butina', 'hierarchical', or 'network'
)
```

**Methods:**
- `analyze_sdf(sdf_path)` - Analyze compounds from SDF file
- `analyze_smiles(smiles_list, names)` - Analyze from SMILES strings
- `visualize(result, output_path, style)` - Generate visualization
- `generate_report(result, output_dir, base_name)` - Generate full report
- `print_summary(result)` - Print text summary

### AnalysisResult

Contains all analysis outputs.

**Attributes:**
- `compounds` - List of Compound objects
- `similarity_matrix` - NxN numpy array of Tanimoto similarities
- `clusters` - List of Cluster objects
- `pathway` - EvolutionPathway with lead compounds
- `summary` - Dictionary of summary statistics

**Methods:**
- `to_dict()` - Convert to serializable dictionary
- `save_json(path)` - Save to JSON file

### EvolutionPathway

The identified optimization pathway.

**Attributes:**
- `initial_leads` - List of initial lead compounds
- `intermediate_leads` - List of intermediate leads
- `final_lead` - The final optimized lead
- `pathway_graph` - NetworkX graph of the pathway

## Algorithm Details

### Lead Identification

The algorithm identifies lead compounds based on:
1. **Cluster centrality**: Leads are structurally central to their clusters
2. **Analog density**: Leads have many close analogs around them
3. **Network position**: Initial leads connect to diverse structures

### Pathway Determination

The optimization pathway is determined by:
1. Scoring compounds for "initial lead" vs "final lead" characteristics
2. Initial leads: high centrality, moderate diversity
3. Final leads: high cluster density, tight analog clusters
4. Connecting leads based on structural similarity

## Testing

```bash
# Run tests
pytest tests/ -v

# With coverage
pytest tests/ --cov=compound_evolution --cov-report=html
```

## Example: Generate Sample Data

```bash
# Generate example kinase inhibitor series
cd examples/
python generate_sample_data.py

# Analyze the sample data
compound-evolution analyze sample_kinase_inhibitors.sdf -o output/
```

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Run tests
5. Submit a pull request

## License

MIT License - see LICENSE file for details.

## Citation

If you use this tool in your research, please cite:

```
Compound Evolution Analyzer: A tool for analyzing lead compound
optimization pathways in drug discovery patents.
```

## Acknowledgments

- RDKit for cheminformatics functionality
- NetworkX for graph algorithms
- Matplotlib for visualizations
