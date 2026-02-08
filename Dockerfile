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
    openpyxl \
    && conda clean -afy

# Copy application code
COPY . .

# Add app to Python path (instead of pip install)
ENV PYTHONPATH=/app

# Expose port
EXPOSE 8080

# Run with gunicorn for production
CMD ["gunicorn", "--bind", "0.0.0.0:8080", "--chdir", "/app", "compound_evolution.web.app:app"]
