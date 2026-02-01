FROM continuumio/miniconda3:latest

WORKDIR /app

# Install RDKit and dependencies via conda
RUN conda install -y -c conda-forge \
    rdkit \
    python=3.10 \
    numpy \
    scipy \
    networkx \
    matplotlib \
    pillow \
    flask \
    gunicorn \
    && conda clean -afy

# Copy application code
COPY . .

# Add app to Python path (instead of pip install)
ENV PYTHONPATH=/app

# Debug: List files to verify they exist
RUN echo "=== Checking template files ===" && \
    ls -la /app/compound_evolution/web/ && \
    ls -la /app/compound_evolution/web/templates/ && \
    ls -la /app/compound_evolution/web/static/

# Expose port
EXPOSE 8080

# Run with gunicorn for production
CMD ["gunicorn", "--bind", "0.0.0.0:8080", "--chdir", "/app", "compound_evolution.web.app:app"]
