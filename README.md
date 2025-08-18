# omics-cancer-explore

Explore TCGA data with this project

## Author
Sivakumar Gowrisankar

---
# TCGA Colorectal Cancer Analysis Pipeline

This project provides an end-to-end bioinformatics pipeline for analyzing colorectal adenocarcinoma (COAD-READ) data from The Cancer Genome Atlas (TCGA). It is designed as a showcase of best practices in data analysis, software engineering, and reproducibility.

The entire pipeline is containerized with Docker and includes unit and integration tests to ensure correctness and reliability.

***

## ## Features

* **Differential Gene Expression (DGE):** Identifies significant gene expression differences between matched tumor and normal tissue samples.
* **Interactive Visualization:** Generates a dynamic and interactive heatmap of the top 100 differentially expressed genes using Plotly.
* **Survival Analysis:** Performs Kaplan-Meier survival analysis to assess the prognostic value of a selected gene's expression on patient survival.
* **Reproducibility:** Packaged with Docker and Conda environments for reproducibility.

***

## ## Project Structure

The project is organized as a standard Python package, separating source code, tests, and configuration for clarity and maintainability.

***

## ## Setup and Usage 🚀

You can run this project using either Docker or a local Conda environment.

### ### Prerequisites

* Git
* Docker (for Option 1)
* Miniconda (for Option 2)

### ### Option 1: Run with Docker

This is the simplest way to run the entire analysis pipeline, as it handles all dependencies within a container.

1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/sivagowrisankar/omics_cancer_explore
    cd omics_cancer_explore
    ```

2.  **Build the Docker Image:**
    Make the build script executable and run it. This will create a Docker image named `tcga_coad_analyzer`.
    ```bash
    chmod +x bin/build_docker.sh
    ./bin/build_docker.sh --outputdir <OutputDir> --survgene <GENE>
    ```

3.  **Run the Analysis Pipeline:**
    Execute the run script. This will start the container, run the `run.py` script, and save all outputs to your local `output/` directory.
    ```bash
    chmod +x bin/run_docker.sh
    ./bin/run_docker.sh --outputdir <Output Directory> --survgene <GENE>
    Example: ./bin/run_docker.sh --outputdir output --survgene ETV4
    ```

### ### Option 2: Run with a Local Conda Environment

1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/sivagowrisankar/omics_cancer_explore
    cd omics_cancer_explore
    ```

2.  **Create and Activate the Conda Environment:**
    ```bash
    conda env create -f docker/omics_cancer_explore/environment.yml
    conda activate omics_cancer_explore
    ```

3.  **Install the Package in Editable Mode:**
    This command links your environment to your source code, making it importable.
    ```bash
    pip install -e .
    ```

4.  **Run the Analysis Pipeline:**
    ```bash
    python omics_cancer_explore/omics_cancer_explore.py --outputdir <Output Directory> --survgene <GENE>
    Example: python omics_cancer_explore/omics_cancer_explore.py --outputdir output --survgene ETV4
    ```

***

## ## Running Tests

The project includes a test suite to verify the functionality of the analysis, and visualization modules.

1.  **Activate Your Conda Environment:**
    ```bash
    conda activate omics_cancer_explore
    ```

2.  **Run `pytest`:**
    Execute `pytest` from the project's root directory. The `-v` flag provides verbose output.
    ```bash
    pytest -v
    ```

***

## ## Expected Outputs

After running the pipeline, the `output/` directory will be populated with the following files:

* `diffexp_heatmap.html`: An interactive Plotly heatmap showing the top 100 differentially expressed genes.
* `{GENE}_surv_curve.png`: A Kaplan-Meier plot illustrating the difference in overall survival between patients with high vs. low expression of the user-given gene.

***

## ## Data Source

The data used in this project is a pre-processed subset of the **TCGA-COADREAD (Colorectal Adenocarcinoma)** cohort, originally sourced from the UCSC Xena dataset for gene expression and the GDC for clinical metadata. The pre-processed CSV files are included in the `data/` directory for immediate use.

***

## ## License

This project is licensed under the MIT License.
