# Start from a lightweight miniconda image
FROM continuumio/miniconda3:24.4.0-0

WORKDIR /app

COPY environment.yml /tmp/environment.yml

RUN conda env create -f /tmp/environment.yml && conda clean -afy

SHELL ["conda", "run", "-n", "omics-cancer-explore", "/bin/bash", "-c"]

COPY . /app

CMD ["bash"]

