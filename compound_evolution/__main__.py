"""
Entry point for running compound_evolution as a module.

Usage:
    python -m compound_evolution analyze input.sdf -o output/
"""

from .cli import main

if __name__ == "__main__":
    main()
