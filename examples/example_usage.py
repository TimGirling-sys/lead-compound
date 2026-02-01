"""
Example usage of the Compound Evolution Analyzer.

This script demonstrates how to use the analyzer programmatically
to analyze compounds and generate visualizations.
"""

import os
import sys

# Add parent directory to path for development
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from compound_evolution import CompoundEvolutionAnalyzer


def example_basic_analysis():
    """
    Basic example: Analyze an SDF file and generate a report.
    """
    print("=" * 60)
    print("Example 1: Basic SDF Analysis")
    print("=" * 60)

    # Initialize analyzer with default settings
    analyzer = CompoundEvolutionAnalyzer(
        similarity_threshold=0.7,  # Tanimoto threshold for analogs
        clustering_method='butina'  # Good for lead identification
    )

    # Path to sample data
    script_dir = os.path.dirname(os.path.abspath(__file__))
    sdf_path = os.path.join(script_dir, "sample_kinase_inhibitors.sdf")

    # Check if sample data exists
    if not os.path.exists(sdf_path):
        print(f"Sample data not found at {sdf_path}")
        print("Run 'python generate_sample_data.py' first to create it.")
        return None

    # Run analysis
    result = analyzer.analyze_sdf(sdf_path)

    # Print summary
    analyzer.print_summary(result)

    return result, analyzer


def example_generate_visualizations(result, analyzer):
    """
    Example: Generate various visualizations from analysis results.
    """
    print("\n" + "=" * 60)
    print("Example 2: Generate Visualizations")
    print("=" * 60)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    # Generate network visualization
    network_path = os.path.join(output_dir, "pathway_network.png")
    analyzer.visualize(
        result,
        network_path,
        style='network',
        title="Kinase Inhibitor Evolution Pathway"
    )
    print(f"Network visualization saved to: {network_path}")

    # Generate linear pathway visualization
    linear_path = os.path.join(output_dir, "pathway_linear.png")
    analyzer.visualize(
        result,
        linear_path,
        style='linear',
        title="Lead Optimization Timeline"
    )
    print(f"Linear visualization saved to: {linear_path}")

    return output_dir


def example_full_report(result, analyzer):
    """
    Example: Generate a complete analysis report with all outputs.
    """
    print("\n" + "=" * 60)
    print("Example 3: Generate Complete Report")
    print("=" * 60)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, "output", "full_report")

    # Generate comprehensive report
    outputs = analyzer.generate_report(
        result,
        output_dir,
        base_name="kinase_inhibitor_evolution"
    )

    print("\nGenerated files:")
    for output_type, path in outputs.items():
        print(f"  {output_type}: {path}")


def example_smiles_analysis():
    """
    Example: Analyze compounds from SMILES strings directly.
    """
    print("\n" + "=" * 60)
    print("Example 4: Analyze SMILES List")
    print("=" * 60)

    # Example SMILES - a small SAR series
    smiles_list = [
        "c1ccc(Nc2ncnc3ccccc23)cc1",           # Parent
        "Fc1ccc(Nc2ncnc3ccccc23)cc1",          # 4-F
        "Fc1cccc(Nc2ncnc3ccccc23)c1",          # 3-F
        "Fc1ccc(Nc2ncnc3ccccc23)cc1F",         # 3,4-diF
        "Clc1ccc(Nc2ncnc3ccccc23)cc1",         # 4-Cl
        "COc1ccc(Nc2ncnc3ccccc23)cc1",         # 4-OMe
    ]

    names = [
        "Parent",
        "4-Fluoro",
        "3-Fluoro",
        "3,4-DiFluoro",
        "4-Chloro",
        "4-Methoxy"
    ]

    analyzer = CompoundEvolutionAnalyzer(
        similarity_threshold=0.65,  # Lower threshold for small series
        clustering_method='network'
    )

    result = analyzer.analyze_smiles(smiles_list, names)
    analyzer.print_summary(result)

    return result


def example_custom_analysis():
    """
    Example: Customize analysis parameters.
    """
    print("\n" + "=" * 60)
    print("Example 5: Custom Analysis Parameters")
    print("=" * 60)

    # Different clustering methods give different results
    print("\nComparing clustering methods:\n")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    sdf_path = os.path.join(script_dir, "sample_kinase_inhibitors.sdf")

    if not os.path.exists(sdf_path):
        print("Sample data not found. Skipping this example.")
        return

    for method in ['butina', 'hierarchical', 'network']:
        analyzer = CompoundEvolutionAnalyzer(
            similarity_threshold=0.7,
            clustering_method=method
        )
        result = analyzer.analyze_sdf(sdf_path)

        print(f"{method.upper()} clustering:")
        print(f"  Clusters found: {len(result.clusters)}")
        print(f"  Leads identified: {len(result.pathway.all_leads)}")
        if result.pathway.final_lead:
            print(f"  Final lead: {result.pathway.final_lead.compound.name}")
        print()


def example_access_raw_data(result):
    """
    Example: Access raw analysis data for custom processing.
    """
    print("\n" + "=" * 60)
    print("Example 6: Access Raw Data")
    print("=" * 60)

    # Access similarity matrix
    print(f"\nSimilarity matrix shape: {result.similarity_matrix.shape}")
    print(f"Average similarity: {result.similarity_matrix.mean():.3f}")

    # Access clusters
    print(f"\nCluster details:")
    for cluster in result.clusters[:3]:  # Show first 3 clusters
        print(f"  Cluster {cluster.id}: {len(cluster.compounds)} compounds")
        print(f"    Centroid: {cluster.centroid.name if cluster.centroid else 'N/A'}")
        print(f"    Avg internal similarity: {cluster.avg_internal_similarity:.3f}")

    # Access pathway leads
    print(f"\nPathway leads:")
    for lead in result.pathway.all_leads:
        print(f"  {lead.lead_type.value.title()}: {lead.compound.name}")
        print(f"    Generation: {lead.generation}")
        print(f"    Analogs: {lead.analog_count}")

    # Export to JSON
    result_dict = result.to_dict()
    print(f"\nResult dictionary keys: {list(result_dict.keys())}")


def main():
    """Run all examples."""
    print("\n" + "#" * 60)
    print("COMPOUND EVOLUTION ANALYZER - EXAMPLES")
    print("#" * 60)

    # Run examples
    output = example_basic_analysis()

    if output is not None:
        result, analyzer = output

        example_generate_visualizations(result, analyzer)
        example_full_report(result, analyzer)
        example_access_raw_data(result)

    example_smiles_analysis()
    example_custom_analysis()

    print("\n" + "#" * 60)
    print("All examples completed!")
    print("#" * 60)


if __name__ == "__main__":
    main()
