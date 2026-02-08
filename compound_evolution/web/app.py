"""
Web application for Compound Evolution Analyzer.
Self-contained version with embedded HTML/CSS/JS - no external templates needed.
"""

import os
import uuid
import shutil
import tempfile
from typing import Optional

from flask import Flask, request, jsonify, send_file, Response
from werkzeug.utils import secure_filename

from ..analyzer import CompoundEvolutionAnalyzer


# Complete HTML with embedded CSS and JavaScript
INDEX_HTML = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Compound Evolution Analyzer</title>
    <style>
        :root {
            --primary-color: #3498db;
            --primary-dark: #2980b9;
            --secondary-color: #2ecc71;
            --danger-color: #e74c3c;
            --text-color: #2c3e50;
            --text-light: #7f8c8d;
            --background: #f5f7fa;
            --card-background: #ffffff;
            --border-color: #e1e8ed;
            --shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
            --radius: 12px;
            --radius-small: 8px;
            --initial-color: #2ecc71;
            --intermediate-color: #3498db;
            --final-color: #e74c3c;
            --analog-color: #95a5a6;
        }
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: var(--background);
            color: var(--text-color);
            line-height: 1.6;
            min-height: 100vh;
        }
        .container { max-width: 1200px; margin: 0 auto; padding: 20px; }
        .header { text-align: center; padding: 40px 20px; }
        .header h1 { font-size: 2.5rem; font-weight: 700; margin-bottom: 10px; }
        .subtitle { font-size: 1.1rem; color: var(--text-light); max-width: 600px; margin: 0 auto; }
        .card {
            background: var(--card-background);
            border-radius: var(--radius);
            padding: 30px;
            margin-bottom: 20px;
            box-shadow: var(--shadow);
        }
        .card h2 { font-size: 1.5rem; margin-bottom: 20px; }
        .upload-area {
            border: 2px dashed var(--border-color);
            border-radius: var(--radius);
            padding: 40px;
            text-align: center;
            transition: all 0.3s ease;
            cursor: pointer;
        }
        .upload-area:hover, .upload-area.drag-over {
            border-color: var(--primary-color);
            background: rgba(52, 152, 219, 0.05);
        }
        .upload-icon { color: var(--text-light); margin-bottom: 15px; }
        .upload-text { font-size: 1.1rem; margin-bottom: 5px; }
        .upload-subtext { color: var(--text-light); margin-bottom: 15px; }
        .file-name { margin-top: 15px; color: var(--primary-color); font-weight: 500; }
        .btn {
            display: inline-flex;
            align-items: center;
            justify-content: center;
            padding: 12px 24px;
            font-size: 1rem;
            font-weight: 500;
            border: none;
            border-radius: var(--radius-small);
            cursor: pointer;
            transition: all 0.2s ease;
            text-decoration: none;
            gap: 8px;
        }
        .btn-primary { background: var(--primary-color); color: white; }
        .btn-primary:hover:not(:disabled) { background: var(--primary-dark); }
        .btn-secondary { background: var(--background); color: var(--text-color); border: 1px solid var(--border-color); }
        .btn-secondary:hover { background: var(--border-color); }
        .btn-outline { background: transparent; color: var(--primary-color); border: 2px solid var(--primary-color); }
        .btn-outline:hover { background: var(--primary-color); color: white; }
        .btn-large { width: 100%; padding: 16px 32px; font-size: 1.1rem; margin-top: 20px; }
        .btn:disabled { opacity: 0.5; cursor: not-allowed; }
        .spinner {
            width: 20px; height: 20px;
            border: 2px solid transparent;
            border-top-color: currentColor;
            border-radius: 50%;
            animation: spin 0.8s linear infinite;
            display: inline-block;
        }
        @keyframes spin { to { transform: rotate(360deg); } }
        .options { margin-top: 30px; padding-top: 20px; border-top: 1px solid var(--border-color); }
        .options h3 { font-size: 1.1rem; margin-bottom: 20px; }
        .option-group { margin-bottom: 20px; }
        .option-group label { display: block; font-weight: 500; margin-bottom: 8px; }
        .option-group select {
            width: 100%; padding: 12px; font-size: 1rem;
            border: 1px solid var(--border-color);
            border-radius: var(--radius-small);
            background: white;
        }
        .option-help { font-size: 0.85rem; color: var(--text-light); margin-top: 5px; }
        .slider-container { display: flex; align-items: center; gap: 15px; }
        .slider-container input[type="range"] {
            flex: 1; height: 6px; -webkit-appearance: none;
            background: var(--border-color); border-radius: 3px; outline: none;
        }
        .slider-container input[type="range"]::-webkit-slider-thumb {
            -webkit-appearance: none; width: 20px; height: 20px;
            background: var(--primary-color); border-radius: 50%; cursor: pointer;
        }
        #threshold-value { min-width: 40px; font-weight: 500; color: var(--primary-color); }
        .progress-container { padding: 20px 0; }
        .progress-bar {
            height: 8px; background: var(--border-color);
            border-radius: 4px; overflow: hidden; margin-bottom: 15px;
        }
        .progress-fill {
            height: 100%;
            background: linear-gradient(90deg, var(--primary-color), var(--secondary-color));
            width: 0%; transition: width 0.3s ease;
        }
        .progress-text { text-align: center; color: var(--text-light); }
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 15px; margin-bottom: 30px;
        }
        .stat-card { background: var(--background); padding: 20px; border-radius: var(--radius-small); text-align: center; }
        .stat-value { font-size: 2rem; font-weight: 700; color: var(--primary-color); }
        .stat-label { font-size: 0.9rem; color: var(--text-light); margin-top: 5px; }
        .tabs { display: flex; gap: 5px; border-bottom: 2px solid var(--border-color); margin-bottom: 20px; overflow-x: auto; }
        .tab {
            padding: 12px 20px; background: none; border: none;
            font-size: 1rem; color: var(--text-light); cursor: pointer;
            position: relative; white-space: nowrap;
        }
        .tab:hover { color: var(--text-color); }
        .tab.active { color: var(--primary-color); font-weight: 500; }
        .tab.active::after {
            content: ''; position: absolute; bottom: -2px; left: 0; right: 0;
            height: 2px; background: var(--primary-color);
        }
        .tab-pane { display: none; }
        .tab-pane.active { display: block; }
        .visualization-container {
            background: white; border: 1px solid var(--border-color);
            border-radius: var(--radius-small); padding: 10px;
            text-align: center; overflow: auto;
        }
        .visualization-container img { max-width: 100%; height: auto; border-radius: var(--radius-small); }
        .visualization-legend { display: flex; justify-content: center; gap: 20px; margin-top: 15px; flex-wrap: wrap; }
        .legend-item { display: flex; align-items: center; gap: 6px; font-size: 0.9rem; color: var(--text-light); }
        .legend-dot { width: 12px; height: 12px; border-radius: 50%; }
        .legend-dot.initial { background: var(--initial-color); }
        .legend-dot.intermediate { background: var(--intermediate-color); }
        .legend-dot.final { background: var(--final-color); }
        .legend-dot.analog { background: var(--analog-color); }
        .pathway-details { display: flex; flex-direction: column; gap: 20px; }
        .pathway-section { border: 1px solid var(--border-color); border-radius: var(--radius-small); overflow: hidden; }
        .pathway-section-header { padding: 15px 20px; font-weight: 600; display: flex; align-items: center; gap: 10px; }
        .pathway-section-header.initial { background: rgba(46, 204, 113, 0.1); color: var(--initial-color); }
        .pathway-section-header.intermediate { background: rgba(52, 152, 219, 0.1); color: var(--intermediate-color); }
        .pathway-section-header.final { background: rgba(231, 76, 60, 0.1); color: var(--final-color); }
        .pathway-section-body { padding: 15px 20px; }
        .lead-info { margin-bottom: 15px; padding-bottom: 15px; border-bottom: 1px solid var(--border-color); }
        .lead-info:last-child { margin-bottom: 0; padding-bottom: 0; border-bottom: none; }
        .lead-name { font-weight: 600; font-size: 1.1rem; margin-bottom: 8px; }
        .lead-smiles {
            font-family: monospace; font-size: 0.85rem; color: var(--text-light);
            word-break: break-all; background: var(--background);
            padding: 8px 12px; border-radius: var(--radius-small); margin-bottom: 10px;
        }
        .lead-meta { display: flex; gap: 20px; font-size: 0.9rem; color: var(--text-light); }
        .pathway-arrow { text-align: center; font-size: 1.5rem; color: var(--text-light); }
        .clusters-table-container { overflow-x: auto; }
        .clusters-table { width: 100%; border-collapse: collapse; }
        .clusters-table th, .clusters-table td { padding: 12px 15px; text-align: left; border-bottom: 1px solid var(--border-color); }
        .clusters-table th { background: var(--background); font-weight: 600; }
        .download-section { margin-top: 30px; padding-top: 20px; border-top: 1px solid var(--border-color); }
        .download-section h3 { margin-bottom: 15px; }
        .download-buttons { display: flex; gap: 10px; flex-wrap: wrap; }
        .new-analysis { margin-top: 30px; text-align: center; }
        .error-card { border-left: 4px solid var(--danger-color); }
        .error-card h2 { color: var(--danger-color); }
        #error-message {
            background: rgba(231, 76, 60, 0.1); padding: 15px;
            border-radius: var(--radius-small); margin-bottom: 20px; color: var(--danger-color);
        }
        .footer { text-align: center; padding: 30px; color: var(--text-light); font-size: 0.9rem; }
        @media (max-width: 768px) {
            .header h1 { font-size: 1.8rem; }
            .card { padding: 20px; }
            .stats-grid { grid-template-columns: repeat(2, 1fr); }
            .download-buttons { flex-direction: column; }
            .download-buttons .btn { width: 100%; }
        }
    </style>
</head>
<body>
    <div class="container">
        <header class="header">
            <h1>Compound Evolution Analyzer</h1>
            <p class="subtitle">Analyze chemical compound evolution in patents to identify lead compounds and optimization pathways</p>
        </header>

        <section id="upload-section" class="card">
            <h2>Upload Compound File</h2>
            <form id="upload-form" enctype="multipart/form-data">
                <div class="upload-area" id="drop-zone">
                    <div class="upload-icon">
                        <svg xmlns="http://www.w3.org/2000/svg" width="64" height="64" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                            <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path>
                            <polyline points="17 8 12 3 7 8"></polyline>
                            <line x1="12" y1="3" x2="12" y2="15"></line>
                        </svg>
                    </div>
                    <p class="upload-text">Drag and drop your file here</p>
                    <p class="upload-subtext">Supported: SDF, Excel (.xlsx, .xls)</p>
                    <label for="file-input" class="btn btn-primary">Choose File</label>
                    <input type="file" id="file-input" name="file" accept=".sdf,.xlsx,.xls" hidden>
                    <p id="file-name" class="file-name"></p>
                </div>

                <div class="options">
                    <h3>Analysis Options</h3>
                    <div class="option-group">
                        <label for="threshold">Similarity Threshold</label>
                        <div class="slider-container">
                            <input type="range" id="threshold" name="threshold" min="0.3" max="0.95" step="0.05" value="0.7">
                            <span id="threshold-value">0.70</span>
                        </div>
                        <p class="option-help">Higher values create stricter analog groupings</p>
                    </div>
                    <div class="option-group">
                        <label for="clustering">Clustering Method</label>
                        <select id="clustering" name="clustering">
                            <option value="butina" selected>Butina (Recommended)</option>
                            <option value="hierarchical">Hierarchical</option>
                            <option value="network">Network-based</option>
                        </select>
                        <p class="option-help">Butina is best for identifying lead series</p>
                    </div>
                </div>

                <button type="submit" class="btn btn-large btn-primary" id="analyze-btn" disabled>
                    <span class="btn-text">Analyze Compounds</span>
                    <span class="btn-loading" hidden><span class="spinner"></span> Analyzing...</span>
                </button>
            </form>
        </section>

        <section id="progress-section" class="card" hidden>
            <h2>Analysis in Progress</h2>
            <div class="progress-container">
                <div class="progress-bar"><div class="progress-fill" id="progress-fill"></div></div>
                <p class="progress-text" id="progress-text">Parsing compounds...</p>
            </div>
        </section>

        <section id="results-section" class="card" hidden>
            <h2>Analysis Results</h2>
            <div id="file-info" style="background: var(--background); padding: 15px; border-radius: var(--radius-small); margin-bottom: 20px;">
                <p style="margin: 0; font-weight: 600;">Analyzed file: <span id="result-filename" style="color: var(--primary-color);"></span></p>
                <p style="margin: 5px 0 0 0; font-size: 0.9rem; color: var(--text-light);">Sample compounds: <span id="sample-compounds"></span></p>
            </div>
            <div class="stats-grid">
                <div class="stat-card"><div class="stat-value" id="stat-compounds">-</div><div class="stat-label">Total Compounds</div></div>
                <div class="stat-card"><div class="stat-value" id="stat-clusters">-</div><div class="stat-label">Clusters</div></div>
                <div class="stat-card"><div class="stat-value" id="stat-leads">-</div><div class="stat-label">Lead Compounds</div></div>
                <div class="stat-card"><div class="stat-value" id="stat-similarity">-</div><div class="stat-label">Avg Similarity</div></div>
                <div class="stat-card"><div class="stat-value" id="stat-molweight">-</div><div class="stat-label">Avg Mol Weight</div></div>
            </div>
            <div class="tabs">
                <button class="tab active" data-tab="network">Network View</button>
                <button class="tab" data-tab="linear">Linear Pathway</button>
                <button class="tab" data-tab="pathway">Pathway Details</button>
                <button class="tab" data-tab="clusters">Clusters</button>
            </div>
            <div class="tab-content">
                <div class="tab-pane active" id="tab-network">
                    <div class="visualization-container"><img id="network-image" src="" alt="Network Visualization"></div>
                    <div class="visualization-legend">
                        <span class="legend-item"><span class="legend-dot initial"></span> Initial Lead</span>
                        <span class="legend-item"><span class="legend-dot intermediate"></span> Intermediate</span>
                        <span class="legend-item"><span class="legend-dot final"></span> Final Lead</span>
                        <span class="legend-item"><span class="legend-dot analog"></span> Analog</span>
                    </div>
                </div>
                <div class="tab-pane" id="tab-linear">
                    <div class="visualization-container"><img id="linear-image" src="" alt="Linear Pathway"></div>
                </div>
                <div class="tab-pane" id="tab-pathway"><div class="pathway-details" id="pathway-details"></div></div>
                <div class="tab-pane" id="tab-clusters">
                    <div class="clusters-table-container">
                        <table class="clusters-table">
                            <thead><tr><th>Cluster</th><th>Size</th><th>Centroid</th><th>Density</th></tr></thead>
                            <tbody id="clusters-tbody"></tbody>
                        </table>
                    </div>
                </div>
            </div>
            <div class="download-section">
                <h3>Download Results</h3>
                <div class="download-buttons">
                    <a id="download-network" class="btn btn-secondary" download>Network PNG</a>
                    <a id="download-linear" class="btn btn-secondary" download>Linear PNG</a>
                    <a id="download-report" class="btn btn-secondary" download>HTML Report</a>
                    <a id="download-data" class="btn btn-secondary" download>JSON Data</a>
                </div>
            </div>
            <div class="new-analysis"><button class="btn btn-outline" id="new-analysis-btn">Start New Analysis</button></div>
        </section>

        <section id="error-section" class="card error-card" hidden>
            <h2>Error</h2>
            <p id="error-message"></p>
            <button class="btn btn-primary" id="try-again-btn">Try Again</button>
        </section>

        <footer class="footer"><p>Compound Evolution Analyzer | Built with RDKit, Flask, and Python</p></footer>
    </div>

    <script>
    document.addEventListener('DOMContentLoaded', function() {
        const uploadSection = document.getElementById('upload-section');
        const progressSection = document.getElementById('progress-section');
        const resultsSection = document.getElementById('results-section');
        const errorSection = document.getElementById('error-section');
        const uploadForm = document.getElementById('upload-form');
        const fileInput = document.getElementById('file-input');
        const dropZone = document.getElementById('drop-zone');
        const fileName = document.getElementById('file-name');
        const analyzeBtn = document.getElementById('analyze-btn');
        const thresholdSlider = document.getElementById('threshold');
        const thresholdValue = document.getElementById('threshold-value');
        const progressFill = document.getElementById('progress-fill');
        const progressText = document.getElementById('progress-text');
        const errorMessage = document.getElementById('error-message');

        let currentAnalysisId = null;
        let selectedFile = null;
        let progressInterval = null;

        dropZone.addEventListener('click', function(e) {
            if (e.target !== fileInput && !e.target.closest('label')) fileInput.click();
        });
        fileInput.addEventListener('change', function(e) { handleFileSelect(e.target.files); });
        dropZone.addEventListener('dragover', function(e) { e.preventDefault(); dropZone.classList.add('drag-over'); });
        dropZone.addEventListener('dragleave', function(e) { e.preventDefault(); dropZone.classList.remove('drag-over'); });
        dropZone.addEventListener('drop', function(e) { e.preventDefault(); dropZone.classList.remove('drag-over'); handleFileSelect(e.dataTransfer.files); });

        function handleFileSelect(files) {
            if (files.length > 0) {
                const file = files[0];
                const name = file.name.toLowerCase();
                const validExts = ['.sdf', '.xlsx', '.xls'];
                if (!validExts.some(ext => name.endsWith(ext))) {
                    alert('Please select an SDF or Excel file (.sdf, .xlsx, .xls)');
                    return;
                }
                selectedFile = file;
                fileName.textContent = file.name;
                analyzeBtn.disabled = false;
            }
        }

        thresholdSlider.addEventListener('input', function() { thresholdValue.textContent = parseFloat(this.value).toFixed(2); });

        document.querySelectorAll('.tab').forEach(function(tab) {
            tab.addEventListener('click', function() {
                document.querySelectorAll('.tab').forEach(function(t) { t.classList.remove('active'); });
                this.classList.add('active');
                document.querySelectorAll('.tab-pane').forEach(function(p) { p.classList.remove('active'); });
                document.getElementById('tab-' + this.dataset.tab).classList.add('active');
            });
        });

        uploadForm.addEventListener('submit', async function(e) {
            e.preventDefault();
            if (!selectedFile) return;
            showSection('progress');
            startProgress();
            const formData = new FormData();
            formData.append('file', selectedFile);
            formData.append('threshold', thresholdSlider.value);
            formData.append('clustering', document.getElementById('clustering').value);
            try {
                const response = await fetch('/upload', { method: 'POST', body: formData });
                const data = await response.json();
                if (!response.ok) throw new Error(data.error || 'Analysis failed');
                currentAnalysisId = data.analysis_id;
                displayResults(data);
                showSection('results');
            } catch (error) {
                errorMessage.textContent = error.message;
                showSection('error');
            }
        });

        document.getElementById('try-again-btn').addEventListener('click', resetAnalysis);
        document.getElementById('new-analysis-btn').addEventListener('click', resetAnalysis);

        function resetAnalysis() {
            if (currentAnalysisId) fetch('/cleanup/' + currentAnalysisId, { method: 'POST' });
            currentAnalysisId = null; selectedFile = null; fileInput.value = '';
            fileName.textContent = ''; analyzeBtn.disabled = true;
            showSection('upload');
        }

        function displayResults(data) {
            document.getElementById('result-filename').textContent = data.filename || 'Unknown';
            document.getElementById('sample-compounds').textContent = (data.summary.sample_compounds || []).join(', ') || 'N/A';
            document.getElementById('stat-compounds').textContent = data.summary.total_compounds;
            document.getElementById('stat-clusters').textContent = data.summary.total_clusters;
            document.getElementById('stat-leads').textContent = data.summary.total_leads;
            document.getElementById('stat-similarity').textContent = data.summary.avg_similarity;
            document.getElementById('stat-molweight').textContent = data.summary.avg_mol_weight || '-';
            document.getElementById('network-image').src = '/image/' + data.analysis_id + '/network?' + Date.now();
            document.getElementById('linear-image').src = '/image/' + data.analysis_id + '/linear?' + Date.now();
            updatePathway(data.pathway);
            updateClusters(data.clusters);
            document.getElementById('download-network').href = '/download/' + data.analysis_id + '/network';
            document.getElementById('download-linear').href = '/download/' + data.analysis_id + '/linear';
            document.getElementById('download-report').href = '/download/' + data.analysis_id + '/report';
            document.getElementById('download-data').href = '/download/' + data.analysis_id + '/data';
        }

        function escapeHtml(text) {
            if (!text) return '';
            var div = document.createElement('div');
            div.textContent = text;
            return div.innerHTML;
        }

        function updatePathway(pathway) {
            var html = '';
            if (pathway.initial_leads && pathway.initial_leads.length > 0) {
                html += '<div class="pathway-section"><div class="pathway-section-header initial"><span class="legend-dot initial"></span>Initial Lead</div><div class="pathway-section-body">';
                pathway.initial_leads.forEach(function(l) { html += '<div class="lead-info"><div class="lead-name">' + escapeHtml(l.name) + '</div><div class="lead-smiles">' + escapeHtml(l.smiles) + '</div><div class="lead-meta"><span>Generation: ' + l.generation + '</span><span>Analogs: ' + l.analog_count + '</span></div></div>'; });
                html += '</div></div>';
            }
            if (pathway.intermediate_leads && pathway.intermediate_leads.length > 0) {
                html += '<div class="pathway-arrow">&#8595;</div><div class="pathway-section"><div class="pathway-section-header intermediate"><span class="legend-dot intermediate"></span>Intermediate</div><div class="pathway-section-body">';
                pathway.intermediate_leads.forEach(function(l) { html += '<div class="lead-info"><div class="lead-name">' + escapeHtml(l.name) + '</div><div class="lead-smiles">' + escapeHtml(l.smiles) + '</div><div class="lead-meta"><span>Generation: ' + l.generation + '</span><span>Analogs: ' + l.analog_count + '</span></div></div>'; });
                html += '</div></div>';
            }
            if (pathway.final_lead) {
                html += '<div class="pathway-arrow">&#8595;</div><div class="pathway-section"><div class="pathway-section-header final"><span class="legend-dot final"></span>Final Lead</div><div class="pathway-section-body">';
                html += '<div class="lead-info"><div class="lead-name">' + escapeHtml(pathway.final_lead.name) + '</div><div class="lead-smiles">' + escapeHtml(pathway.final_lead.smiles) + '</div><div class="lead-meta"><span>Generation: ' + pathway.final_lead.generation + '</span><span>Analogs: ' + pathway.final_lead.analog_count + '</span></div></div>';
                html += '</div></div>';
            }
            document.getElementById('pathway-details').innerHTML = html;
        }

        function updateClusters(clusters) {
            document.getElementById('clusters-tbody').innerHTML = clusters.map(function(c) {
                return '<tr><td>Cluster ' + c.id + '</td><td>' + c.size + '</td><td>' + escapeHtml(c.centroid || 'N/A') + '</td><td>' + c.density + '</td></tr>';
            }).join('');
        }

        function startProgress() {
            var progress = 0;
            var messages = ['Parsing compounds...', 'Calculating fingerprints...', 'Computing similarity...', 'Clustering compounds...', 'Identifying leads...', 'Generating visualizations...'];
            progressInterval = setInterval(function() {
                progress += Math.random() * 15;
                if (progress > 95) progress = 95;
                progressFill.style.width = progress + '%';
                progressText.textContent = messages[Math.min(Math.floor(progress / 16), messages.length - 1)];
            }, 500);
        }

        function showSection(name) {
            if (progressInterval) { clearInterval(progressInterval); progressInterval = null; }
            uploadSection.hidden = name !== 'upload';
            progressSection.hidden = name !== 'progress';
            resultsSection.hidden = name !== 'results';
            errorSection.hidden = name !== 'error';
            analyzeBtn.querySelector('.btn-text').hidden = name === 'progress';
            analyzeBtn.querySelector('.btn-loading').hidden = name !== 'progress';
        }
    });
    </script>
</body>
</html>'''


def create_app(config: Optional[dict] = None) -> Flask:
    """Create and configure the Flask application."""
    app = Flask(__name__)

    app.config.update(
        SECRET_KEY=os.environ.get('SECRET_KEY', 'dev-key-change-in-production'),
        MAX_CONTENT_LENGTH=50 * 1024 * 1024,
        UPLOAD_FOLDER=tempfile.mkdtemp(prefix='compound_evolution_'),
        ALLOWED_EXTENSIONS={'sdf', 'xlsx', 'xls'},
    )

    if config:
        app.config.update(config)

    os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
    register_routes(app)
    return app


def allowed_file(filename: str, allowed_extensions: set) -> bool:
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in allowed_extensions


def register_routes(app: Flask) -> None:

    @app.route('/')
    def index():
        return Response(INDEX_HTML, mimetype='text/html')

    @app.route('/upload', methods=['POST'])
    def upload_file():
        if 'file' not in request.files:
            return jsonify({'error': 'No file provided'}), 400
        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'No file selected'}), 400
        if not allowed_file(file.filename, app.config['ALLOWED_EXTENSIONS']):
            return jsonify({'error': 'Invalid file type. Please upload an SDF or Excel file (.sdf, .xlsx, .xls).'}), 400

        analysis_id = str(uuid.uuid4())
        analysis_dir = os.path.join(app.config['UPLOAD_FOLDER'], analysis_id)
        os.makedirs(analysis_dir, exist_ok=True)

        try:
            filename = secure_filename(file.filename)
            sdf_path = os.path.join(analysis_dir, filename)
            file.save(sdf_path)

            similarity_threshold = float(request.form.get('threshold', 0.7))
            clustering_method = request.form.get('clustering', 'butina')
            similarity_threshold = max(0.3, min(0.95, similarity_threshold))
            if clustering_method not in ['butina', 'hierarchical', 'network']:
                clustering_method = 'butina'

            analyzer = CompoundEvolutionAnalyzer(
                similarity_threshold=similarity_threshold,
                clustering_method=clustering_method
            )
            result = analyzer.analyze_file(sdf_path)

            network_path = os.path.join(analysis_dir, 'network.png')
            linear_path = os.path.join(analysis_dir, 'linear.png')
            analyzer.visualize(result, network_path, style='network')
            analyzer.visualize(result, linear_path, style='linear')

            report_path = os.path.join(analysis_dir, 'report.html')
            analyzer.visualizer.create_summary_report(result.pathway, result.compounds, result.clusters, report_path)

            json_path = os.path.join(analysis_dir, 'data.json')
            result.save_json(json_path)

            # Get sample compound names to prove file was analyzed
            sample_compounds = [c.name for c in result.compounds[:5]]

            return jsonify({
                'success': True,
                'analysis_id': analysis_id,
                'filename': filename,
                'summary': {
                    'total_compounds': len(result.compounds),
                    'total_clusters': len(result.clusters),
                    'total_leads': len(result.pathway.all_leads),
                    'avg_similarity': round(result.summary['avg_pairwise_similarity'], 3),
                    'avg_mol_weight': round(result.summary.get('avg_mol_weight', 0), 1),
                    'sample_compounds': sample_compounds,
                },
                'pathway': {
                    'initial_leads': [{'name': l.compound.name, 'smiles': l.compound.smiles, 'analog_count': l.analog_count, 'generation': l.generation} for l in result.pathway.initial_leads],
                    'intermediate_leads': [{'name': l.compound.name, 'smiles': l.compound.smiles, 'analog_count': l.analog_count, 'generation': l.generation} for l in result.pathway.intermediate_leads],
                    'final_lead': {'name': result.pathway.final_lead.compound.name, 'smiles': result.pathway.final_lead.compound.smiles, 'analog_count': result.pathway.final_lead.analog_count, 'generation': result.pathway.final_lead.generation} if result.pathway.final_lead else None
                },
                'clusters': [{'id': c.id, 'size': len(c.compounds), 'centroid': c.centroid.name if c.centroid else None, 'density': round(c.density, 3)} for c in result.clusters[:10]]
            })
        except Exception as e:
            if os.path.exists(analysis_dir):
                shutil.rmtree(analysis_dir)
            return jsonify({'error': str(e)}), 500

    @app.route('/image/<analysis_id>/<image_type>')
    def get_image(analysis_id: str, image_type: str):
        if image_type not in ['network', 'linear']:
            return jsonify({'error': 'Invalid image type'}), 400
        try:
            uuid.UUID(analysis_id)
        except ValueError:
            return jsonify({'error': 'Invalid analysis ID'}), 400
        image_path = os.path.join(app.config['UPLOAD_FOLDER'], analysis_id, f'{image_type}.png')
        if not os.path.exists(image_path):
            return jsonify({'error': 'Image not found'}), 404
        return send_file(image_path, mimetype='image/png')

    @app.route('/download/<analysis_id>/<file_type>')
    def download_file(analysis_id: str, file_type: str):
        file_map = {'network': ('network.png', 'image/png'), 'linear': ('linear.png', 'image/png'), 'report': ('report.html', 'text/html'), 'data': ('data.json', 'application/json')}
        if file_type not in file_map:
            return jsonify({'error': 'Invalid file type'}), 400
        try:
            uuid.UUID(analysis_id)
        except ValueError:
            return jsonify({'error': 'Invalid analysis ID'}), 400
        filename, mimetype = file_map[file_type]
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], analysis_id, filename)
        if not os.path.exists(file_path):
            return jsonify({'error': 'File not found'}), 404
        return send_file(file_path, mimetype=mimetype, as_attachment=True, download_name=f'compound_evolution_{file_type}.{filename.split(".")[-1]}')

    @app.route('/cleanup/<analysis_id>', methods=['POST'])
    def cleanup(analysis_id: str):
        try:
            uuid.UUID(analysis_id)
        except ValueError:
            return jsonify({'error': 'Invalid analysis ID'}), 400
        analysis_dir = os.path.join(app.config['UPLOAD_FOLDER'], analysis_id)
        if os.path.exists(analysis_dir):
            shutil.rmtree(analysis_dir)
        return jsonify({'success': True})


app = create_app()


def run_server(host: str = '127.0.0.1', port: int = 5000, debug: bool = False):
    print(f"\n{'='*60}")
    print("Compound Evolution Analyzer - Web Interface")
    print(f"{'='*60}")
    print(f"\nStarting server at http://{host}:{port}")
    print("Press Ctrl+C to stop\n")
    app.run(host=host, port=port, debug=debug)


if __name__ == '__main__':
    run_server(debug=True)
