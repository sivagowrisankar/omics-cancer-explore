# Start from a lightweight miniconda image
FROM continuumio/miniconda3:24.4.0-0

WORKDIR /app

COPY docker/omics_cancer_explore/environment.yml .

RUN conda env create -f environment.yml && conda clean -afy

COPY omics_cancer_explore/ ./omics_cancer_explore/
COPY data/ ./data/

SHELL ["conda", "run", "-n", "omics_cancer_explore", "/bin/bash", "-c"]

CMD ["python", "omics_cancer_explore/omics_cancer_explore.py"]

