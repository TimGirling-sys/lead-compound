"""
Web application for Compound Evolution Analyzer.

Provides a user-friendly web interface for uploading SDF files
and visualizing compound evolution pathways.
"""

import os
import uuid
import shutil
import tempfile
from pathlib import Path
from typing import Optional
import base64
import io

from flask import (
    Flask, render_template, request, jsonify,
    send_file, session, redirect, url_for
)
from werkzeug.utils import secure_filename

from ..analyzer import CompoundEvolutionAnalyzer, AnalysisResult
from ..visualization import EvolutionVisualizer


def create_app(config: Optional[dict] = None) -> Flask:
    """
    Create and configure the Flask application.

    Args:
        config: Optional configuration dictionary

    Returns:
        Configured Flask application
    """
    # Get the directory where this file is located
    base_dir = os.path.dirname(os.path.abspath(__file__))
    template_dir = os.path.join(base_dir, 'templates')
    static_dir = os.path.join(base_dir, 'static')

    app = Flask(
        __name__,
        template_folder=template_dir,
        static_folder=static_dir
    )

    # Default configuration
    app.config.update(
        SECRET_KEY=os.environ.get('SECRET_KEY', 'dev-key-change-in-production'),
        MAX_CONTENT_LENGTH=50 * 1024 * 1024,  # 50MB max file size
        UPLOAD_FOLDER=tempfile.mkdtemp(prefix='compound_evolution_'),
        ALLOWED_EXTENSIONS={'sdf'},
    )

    # Override with provided config
    if config:
        app.config.update(config)

    # Ensure upload folder exists
    os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

    # Register routes
    register_routes(app)

    return app


def allowed_file(filename: str, allowed_extensions: set) -> bool:
    """Check if file extension is allowed."""
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in allowed_extensions


def register_routes(app: Flask) -> None:
    """Register all application routes."""

    @app.route('/')
    def index():
        """Main page with upload form."""
        return render_template('index.html')

    @app.route('/upload', methods=['POST'])
    def upload_file():
        """Handle SDF file upload and analysis."""
        if 'file' not in request.files:
            return jsonify({'error': 'No file provided'}), 400

        file = request.files['file']

        if file.filename == '':
            return jsonify({'error': 'No file selected'}), 400

        if not allowed_file(file.filename, app.config['ALLOWED_EXTENSIONS']):
            return jsonify({'error': 'Invalid file type. Please upload an SDF file.'}), 400

        # Generate unique session ID for this analysis
        analysis_id = str(uuid.uuid4())
        analysis_dir = os.path.join(app.config['UPLOAD_FOLDER'], analysis_id)
        os.makedirs(analysis_dir, exist_ok=True)

        try:
            # Save uploaded file
            filename = secure_filename(file.filename)
            sdf_path = os.path.join(analysis_dir, filename)
            file.save(sdf_path)

            # Get analysis parameters from form
            similarity_threshold = float(request.form.get('threshold', 0.7))
            clustering_method = request.form.get('clustering', 'butina')

            # Validate parameters
            similarity_threshold = max(0.3, min(0.95, similarity_threshold))
            if clustering_method not in ['butina', 'hierarchical', 'network']:
                clustering_method = 'butina'

            # Run analysis
            analyzer = CompoundEvolutionAnalyzer(
                similarity_threshold=similarity_threshold,
                clustering_method=clustering_method
            )

            result = analyzer.analyze_sdf(sdf_path)

            # Generate visualizations
            network_path = os.path.join(analysis_dir, 'network.png')
            linear_path = os.path.join(analysis_dir, 'linear.png')

            analyzer.visualize(result, network_path, style='network')
            analyzer.visualize(result, linear_path, style='linear')

            # Generate HTML report
            report_path = os.path.join(analysis_dir, 'report.html')
            analyzer.visualizer.create_summary_report(
                result.pathway,
                result.compounds,
                result.clusters,
                report_path
            )

            # Save JSON data
            json_path = os.path.join(analysis_dir, 'data.json')
            result.save_json(json_path)

            # Prepare response data
            response_data = {
                'success': True,
                'analysis_id': analysis_id,
                'summary': {
                    'total_compounds': len(result.compounds),
                    'total_clusters': len(result.clusters),
                    'total_leads': len(result.pathway.all_leads),
                    'avg_similarity': round(result.summary['avg_pairwise_similarity'], 3),
                },
                'pathway': {
                    'initial_leads': [
                        {
                            'name': l.compound.name,
                            'smiles': l.compound.smiles,
                            'analog_count': l.analog_count,
                            'generation': l.generation
                        }
                        for l in result.pathway.initial_leads
                    ],
                    'intermediate_leads': [
                        {
                            'name': l.compound.name,
                            'smiles': l.compound.smiles,
                            'analog_count': l.analog_count,
                            'generation': l.generation
                        }
                        for l in result.pathway.intermediate_leads
                    ],
                    'final_lead': {
                        'name': result.pathway.final_lead.compound.name,
                        'smiles': result.pathway.final_lead.compound.smiles,
                        'analog_count': result.pathway.final_lead.analog_count,
                        'generation': result.pathway.final_lead.generation
                    } if result.pathway.final_lead else None
                },
                'clusters': [
                    {
                        'id': c.id,
                        'size': len(c.compounds),
                        'centroid': c.centroid.name if c.centroid else None,
                        'density': round(c.density, 3)
                    }
                    for c in result.clusters[:10]  # Top 10 clusters
                ]
            }

            return jsonify(response_data)

        except Exception as e:
            # Clean up on error
            if os.path.exists(analysis_dir):
                shutil.rmtree(analysis_dir)
            return jsonify({'error': str(e)}), 500

    @app.route('/image/<analysis_id>/<image_type>')
    def get_image(analysis_id: str, image_type: str):
        """Serve generated visualization images."""
        if image_type not in ['network', 'linear']:
            return jsonify({'error': 'Invalid image type'}), 400

        # Validate analysis_id format (UUID)
        try:
            uuid.UUID(analysis_id)
        except ValueError:
            return jsonify({'error': 'Invalid analysis ID'}), 400

        image_path = os.path.join(
            app.config['UPLOAD_FOLDER'],
            analysis_id,
            f'{image_type}.png'
        )

        if not os.path.exists(image_path):
            return jsonify({'error': 'Image not found'}), 404

        return send_file(image_path, mimetype='image/png')

    @app.route('/report/<analysis_id>')
    def get_report(analysis_id: str):
        """Serve the HTML report."""
        try:
            uuid.UUID(analysis_id)
        except ValueError:
            return jsonify({'error': 'Invalid analysis ID'}), 400

        report_path = os.path.join(
            app.config['UPLOAD_FOLDER'],
            analysis_id,
            'report.html'
        )

        if not os.path.exists(report_path):
            return jsonify({'error': 'Report not found'}), 404

        return send_file(report_path, mimetype='text/html')

    @app.route('/download/<analysis_id>/<file_type>')
    def download_file(analysis_id: str, file_type: str):
        """Download generated files."""
        file_map = {
            'network': ('network.png', 'image/png'),
            'linear': ('linear.png', 'image/png'),
            'report': ('report.html', 'text/html'),
            'data': ('data.json', 'application/json'),
        }

        if file_type not in file_map:
            return jsonify({'error': 'Invalid file type'}), 400

        try:
            uuid.UUID(analysis_id)
        except ValueError:
            return jsonify({'error': 'Invalid analysis ID'}), 400

        filename, mimetype = file_map[file_type]
        file_path = os.path.join(
            app.config['UPLOAD_FOLDER'],
            analysis_id,
            filename
        )

        if not os.path.exists(file_path):
            return jsonify({'error': 'File not found'}), 404

        return send_file(
            file_path,
            mimetype=mimetype,
            as_attachment=True,
            download_name=f'compound_evolution_{file_type}.{filename.split(".")[-1]}'
        )

    @app.route('/cleanup/<analysis_id>', methods=['POST'])
    def cleanup(analysis_id: str):
        """Clean up analysis files."""
        try:
            uuid.UUID(analysis_id)
        except ValueError:
            return jsonify({'error': 'Invalid analysis ID'}), 400

        analysis_dir = os.path.join(app.config['UPLOAD_FOLDER'], analysis_id)
        if os.path.exists(analysis_dir):
            shutil.rmtree(analysis_dir)

        return jsonify({'success': True})

    @app.errorhandler(413)
    def too_large(e):
        """Handle file too large error."""
        return jsonify({
            'error': 'File too large. Maximum size is 50MB.'
        }), 413

    @app.errorhandler(500)
    def server_error(e):
        """Handle internal server errors."""
        return jsonify({
            'error': 'An internal error occurred. Please try again.'
        }), 500


# Create default app instance
app = create_app()


def run_server(host: str = '127.0.0.1', port: int = 5000, debug: bool = False):
    """
    Run the web server.

    Args:
        host: Host address to bind to
        port: Port number
        debug: Enable debug mode
    """
    print(f"\n{'='*60}")
    print("Compound Evolution Analyzer - Web Interface")
    print(f"{'='*60}")
    print(f"\nStarting server at http://{host}:{port}")
    print("Press Ctrl+C to stop\n")

    app.run(host=host, port=port, debug=debug)


if __name__ == '__main__':
    run_server(debug=True)
