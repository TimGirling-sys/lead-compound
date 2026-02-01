"""
Web interface for Compound Evolution Analyzer.

Provides a Flask-based web application for uploading SDF files
and visualizing compound evolution pathways.
"""

from .app import create_app, run_server, app

__all__ = ['create_app', 'run_server', 'app']
