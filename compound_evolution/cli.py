"""
Command-line interface for compound evolution analysis.

Provides a user-friendly CLI for analyzing patent compounds and
generating evolution pathway visualizations.
"""

import argparse
import sys
from pathlib import Path

from .analyzer import CompoundEvolutionAnalyzer


def create_parser() -> argparse.ArgumentParser:
    """Create the argument parser for the CLI."""
    parser = argparse.ArgumentParser(
        prog='compound-evolution',
        description="""
Analyze chemical compound evolution in patents to identify lead compounds
and their optimization pathways.

This tool parses SDF files containing patent compounds, calculates
structural similarities, identifies clusters of analogs, and determines
the progression from initial lead to final optimized lead compound.
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis with default settings
  compound-evolution analyze patent_compounds.sdf -o results/

  # Analysis with custom similarity threshold
  compound-evolution analyze compounds.sdf -o output/ --threshold 0.8

  # Use different clustering method
  compound-evolution analyze compounds.sdf -o output/ --clustering network

  # Generate only network visualization
  compound-evolution visualize results/evolution_analysis_data.json \\
      -o pathway.png --style network
        """
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Analyze command
    analyze_parser = subparsers.add_parser(
        'analyze',
        help='Analyze compounds from an SDF file'
    )
    analyze_parser.add_argument(
        'input',
        type=str,
        help='Path to input SDF file'
    )
    analyze_parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='Output directory for results'
    )
    analyze_parser.add_argument(
        '--threshold',
        type=float,
        default=0.7,
        help='Tanimoto similarity threshold (default: 0.7)'
    )
    analyze_parser.add_argument(
        '--clustering',
        type=str,
        choices=['hierarchical', 'network', 'butina'],
        default='butina',
        help='Clustering method (default: butina)'
    )
    analyze_parser.add_argument(
        '--fp-radius',
        type=int,
        default=2,
        help='Morgan fingerprint radius (default: 2 = ECFP4)'
    )
    analyze_parser.add_argument(
        '--fp-bits',
        type=int,
        default=2048,
        help='Fingerprint bit vector length (default: 2048)'
    )
    analyze_parser.add_argument(
        '--no-visualization',
        action='store_true',
        help='Skip visualization generation'
    )
    analyze_parser.add_argument(
        '--name',
        type=str,
        default='evolution_analysis',
        help='Base name for output files (default: evolution_analysis)'
    )
    analyze_parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print detailed progress information'
    )

    # Visualize command (for regenerating visualizations)
    viz_parser = subparsers.add_parser(
        'visualize',
        help='Generate visualization from previous analysis'
    )
    viz_parser.add_argument(
        'input',
        type=str,
        help='Path to analysis JSON file'
    )
    viz_parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='Output path for visualization'
    )
    viz_parser.add_argument(
        '--style',
        type=str,
        choices=['network', 'linear'],
        default='network',
        help='Visualization style (default: network)'
    )
    viz_parser.add_argument(
        '--dpi',
        type=int,
        default=150,
        help='Output resolution (default: 150)'
    )
    viz_parser.add_argument(
        '--figsize',
        type=str,
        default='16,12',
        help='Figure size as width,height (default: 16,12)'
    )

    return parser


def cmd_analyze(args: argparse.Namespace) -> int:
    """Execute the analyze command."""
    # Validate input file
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        return 1

    if not input_path.suffix.lower() == '.sdf':
        print(f"Warning: Input file does not have .sdf extension")

    # Create output directory
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)

    try:
        # Initialize analyzer
        analyzer = CompoundEvolutionAnalyzer(
            similarity_threshold=args.threshold,
            fingerprint_radius=args.fp_radius,
            fingerprint_bits=args.fp_bits,
            clustering_method=args.clustering
        )

        # Run analysis
        print(f"\nAnalyzing {args.input}...")
        print(f"  Similarity threshold: {args.threshold}")
        print(f"  Clustering method: {args.clustering}")
        print(f"  Fingerprint: ECFP{args.fp_radius * 2} ({args.fp_bits} bits)\n")

        result = analyzer.analyze_sdf(str(input_path))

        # Print summary
        analyzer.print_summary(result)

        # Generate outputs
        if not args.no_visualization:
            print("\nGenerating outputs...")
            outputs = analyzer.generate_report(
                result,
                str(output_path),
                base_name=args.name
            )
        else:
            # Just save JSON
            json_path = output_path / f"{args.name}_data.json"
            result.save_json(str(json_path))
            print(f"\nSaved analysis data to {json_path}")

        print("\nAnalysis complete!")
        return 0

    except Exception as e:
        print(f"Error during analysis: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


def cmd_visualize(args: argparse.Namespace) -> int:
    """Execute the visualize command."""
    print("Note: Visualization from JSON not yet implemented.")
    print("Please use the 'analyze' command to generate visualizations.")
    return 1


def main() -> int:
    """Main entry point for the CLI."""
    parser = create_parser()
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return 0

    if args.command == 'analyze':
        return cmd_analyze(args)
    elif args.command == 'visualize':
        return cmd_visualize(args)
    else:
        parser.print_help()
        return 1


if __name__ == '__main__':
    sys.exit(main())
