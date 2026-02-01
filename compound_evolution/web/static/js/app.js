/**
 * Compound Evolution Analyzer - Frontend JavaScript
 */

document.addEventListener('DOMContentLoaded', function() {
    // Elements
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
    const tryAgainBtn = document.getElementById('try-again-btn');
    const newAnalysisBtn = document.getElementById('new-analysis-btn');

    // State
    let currentAnalysisId = null;
    let selectedFile = null;

    // Initialize
    init();

    function init() {
        setupFileUpload();
        setupSlider();
        setupTabs();
        setupFormSubmission();
        setupResetButtons();
    }

    // File Upload Setup
    function setupFileUpload() {
        // Click to upload
        dropZone.addEventListener('click', function(e) {
            if (e.target !== fileInput && !e.target.closest('label')) {
                fileInput.click();
            }
        });

        // File selected via input
        fileInput.addEventListener('change', function(e) {
            handleFileSelect(e.target.files);
        });

        // Drag and drop
        dropZone.addEventListener('dragover', function(e) {
            e.preventDefault();
            dropZone.classList.add('drag-over');
        });

        dropZone.addEventListener('dragleave', function(e) {
            e.preventDefault();
            dropZone.classList.remove('drag-over');
        });

        dropZone.addEventListener('drop', function(e) {
            e.preventDefault();
            dropZone.classList.remove('drag-over');
            handleFileSelect(e.dataTransfer.files);
        });
    }

    function handleFileSelect(files) {
        if (files.length > 0) {
            const file = files[0];

            // Validate file type
            if (!file.name.toLowerCase().endsWith('.sdf')) {
                showError('Please select an SDF file.');
                return;
            }

            selectedFile = file;
            fileName.textContent = file.name;
            analyzeBtn.disabled = false;
        }
    }

    // Slider Setup
    function setupSlider() {
        thresholdSlider.addEventListener('input', function() {
            thresholdValue.textContent = parseFloat(this.value).toFixed(2);
        });
    }

    // Tabs Setup
    function setupTabs() {
        const tabs = document.querySelectorAll('.tab');
        tabs.forEach(tab => {
            tab.addEventListener('click', function() {
                const targetId = this.dataset.tab;

                // Update active tab
                tabs.forEach(t => t.classList.remove('active'));
                this.classList.add('active');

                // Update active pane
                document.querySelectorAll('.tab-pane').forEach(pane => {
                    pane.classList.remove('active');
                });
                document.getElementById('tab-' + targetId).classList.add('active');
            });
        });
    }

    // Form Submission
    function setupFormSubmission() {
        uploadForm.addEventListener('submit', async function(e) {
            e.preventDefault();

            if (!selectedFile) {
                showError('Please select a file first.');
                return;
            }

            // Show progress
            showSection('progress');
            startProgressAnimation();

            // Prepare form data
            const formData = new FormData();
            formData.append('file', selectedFile);
            formData.append('threshold', thresholdSlider.value);
            formData.append('clustering', document.getElementById('clustering').value);

            try {
                const response = await fetch('/upload', {
                    method: 'POST',
                    body: formData
                });

                const data = await response.json();

                if (!response.ok) {
                    throw new Error(data.error || 'Analysis failed');
                }

                // Store analysis ID
                currentAnalysisId = data.analysis_id;

                // Show results
                displayResults(data);
                showSection('results');

            } catch (error) {
                showError(error.message);
                showSection('error');
            }
        });
    }

    // Reset Buttons
    function setupResetButtons() {
        tryAgainBtn.addEventListener('click', resetAnalysis);
        newAnalysisBtn.addEventListener('click', resetAnalysis);
    }

    function resetAnalysis() {
        // Clean up previous analysis
        if (currentAnalysisId) {
            fetch(`/cleanup/${currentAnalysisId}`, { method: 'POST' })
                .catch(() => {}); // Ignore errors
        }

        // Reset state
        currentAnalysisId = null;
        selectedFile = null;
        fileInput.value = '';
        fileName.textContent = '';
        analyzeBtn.disabled = true;

        // Show upload section
        showSection('upload');
    }

    // Display Results
    function displayResults(data) {
        // Update stats
        document.getElementById('stat-compounds').textContent = data.summary.total_compounds;
        document.getElementById('stat-clusters').textContent = data.summary.total_clusters;
        document.getElementById('stat-leads').textContent = data.summary.total_leads;
        document.getElementById('stat-similarity').textContent = data.summary.avg_similarity;

        // Update images
        document.getElementById('network-image').src = `/image/${data.analysis_id}/network?t=${Date.now()}`;
        document.getElementById('linear-image').src = `/image/${data.analysis_id}/linear?t=${Date.now()}`;

        // Update pathway details
        updatePathwayDetails(data.pathway);

        // Update clusters table
        updateClustersTable(data.clusters);

        // Update download links
        document.getElementById('download-network').href = `/download/${data.analysis_id}/network`;
        document.getElementById('download-linear').href = `/download/${data.analysis_id}/linear`;
        document.getElementById('download-report').href = `/download/${data.analysis_id}/report`;
        document.getElementById('download-data').href = `/download/${data.analysis_id}/data`;
    }

    function updatePathwayDetails(pathway) {
        const container = document.getElementById('pathway-details');
        let html = '';

        // Initial leads
        if (pathway.initial_leads && pathway.initial_leads.length > 0) {
            html += `
                <div class="pathway-section">
                    <div class="pathway-section-header initial">
                        <span class="legend-dot initial"></span>
                        Initial Lead${pathway.initial_leads.length > 1 ? 's' : ''}
                    </div>
                    <div class="pathway-section-body">
                        ${pathway.initial_leads.map(lead => createLeadHtml(lead)).join('')}
                    </div>
                </div>
            `;
        }

        // Arrow
        if (pathway.intermediate_leads && pathway.intermediate_leads.length > 0) {
            html += '<div class="pathway-arrow">&#8595;</div>';
        }

        // Intermediate leads
        if (pathway.intermediate_leads && pathway.intermediate_leads.length > 0) {
            html += `
                <div class="pathway-section">
                    <div class="pathway-section-header intermediate">
                        <span class="legend-dot intermediate"></span>
                        Intermediate Lead${pathway.intermediate_leads.length > 1 ? 's' : ''}
                    </div>
                    <div class="pathway-section-body">
                        ${pathway.intermediate_leads.map(lead => createLeadHtml(lead)).join('')}
                    </div>
                </div>
            `;
        }

        // Arrow
        if (pathway.final_lead) {
            html += '<div class="pathway-arrow">&#8595;</div>';
        }

        // Final lead
        if (pathway.final_lead) {
            html += `
                <div class="pathway-section">
                    <div class="pathway-section-header final">
                        <span class="legend-dot final"></span>
                        Final Lead
                    </div>
                    <div class="pathway-section-body">
                        ${createLeadHtml(pathway.final_lead)}
                    </div>
                </div>
            `;
        }

        container.innerHTML = html;
    }

    function createLeadHtml(lead) {
        return `
            <div class="lead-info">
                <div class="lead-name">${escapeHtml(lead.name)}</div>
                <div class="lead-smiles">${escapeHtml(lead.smiles)}</div>
                <div class="lead-meta">
                    <span>Generation: ${lead.generation}</span>
                    <span>Analogs: ${lead.analog_count}</span>
                </div>
            </div>
        `;
    }

    function updateClustersTable(clusters) {
        const tbody = document.getElementById('clusters-tbody');
        tbody.innerHTML = clusters.map(cluster => `
            <tr>
                <td>Cluster ${cluster.id}</td>
                <td>${cluster.size}</td>
                <td>${escapeHtml(cluster.centroid || 'N/A')}</td>
                <td>${cluster.density}</td>
            </tr>
        `).join('');
    }

    // Progress Animation
    function startProgressAnimation() {
        let progress = 0;
        const messages = [
            'Parsing compounds...',
            'Calculating fingerprints...',
            'Computing similarity matrix...',
            'Clustering compounds...',
            'Identifying lead compounds...',
            'Analyzing pathway...',
            'Generating visualizations...',
            'Finalizing results...'
        ];

        const interval = setInterval(() => {
            progress += Math.random() * 15;
            if (progress > 95) progress = 95;

            progressFill.style.width = progress + '%';

            const messageIndex = Math.min(
                Math.floor(progress / 12.5),
                messages.length - 1
            );
            progressText.textContent = messages[messageIndex];
        }, 500);

        // Store interval ID for cleanup
        progressFill.dataset.interval = interval;
    }

    function stopProgressAnimation() {
        const interval = progressFill.dataset.interval;
        if (interval) {
            clearInterval(parseInt(interval));
            progressFill.style.width = '100%';
        }
    }

    // Section Management
    function showSection(sectionName) {
        stopProgressAnimation();

        uploadSection.hidden = sectionName !== 'upload';
        progressSection.hidden = sectionName !== 'progress';
        resultsSection.hidden = sectionName !== 'results';
        errorSection.hidden = sectionName !== 'error';

        // Reset button state
        if (sectionName !== 'progress') {
            analyzeBtn.querySelector('.btn-text').hidden = false;
            analyzeBtn.querySelector('.btn-loading').hidden = true;
        } else {
            analyzeBtn.querySelector('.btn-text').hidden = true;
            analyzeBtn.querySelector('.btn-loading').hidden = false;
        }
    }

    function showError(message) {
        errorMessage.textContent = message;
    }

    // Utility
    function escapeHtml(text) {
        if (!text) return '';
        const div = document.createElement('div');
        div.textContent = text;
        return div.innerHTML;
    }
});
