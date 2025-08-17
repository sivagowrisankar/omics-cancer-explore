#!/bin/bash
# A script to run the analysis inside the Docker container, with argument handling.

set -e

IMAGE_NAME="omics_cancer_analyzer"
TAG="latest"

# --- Argument Parsing ---
# Default host output directory
HOST_OUTPUT_DIR="output"
# Array to hold arguments
CONTAINER_ARGS=()

# Separate --outputdir for volume mapping from other args
while [[ $# -gt 0 ]]; do
  case $1 in
    --outputdir)
      HOST_OUTPUT_DIR="$2"
      # We still pass the argument to the container so it knows where to write
      CONTAINER_ARGS+=("$1" "$2")
      shift # past argument
      shift # past value
      ;;
    *)
      CONTAINER_ARGS+=("$1") # save other args
      shift # past argument
      ;;
  esac
done

# --- Execution ---
# Get the absolute path of the project directory
PROJECT_DIR=$(pwd)
HOST_OUTPUT_PATH="${PROJECT_DIR}/${HOST_OUTPUT_DIR}"

# Create the output directory on the host if it doesn't exist
mkdir -p "${HOST_OUTPUT_PATH}"

echo "Running analysis in Docker container..."
echo "Host output directory: ${HOST_OUTPUT_PATH}"
echo "Passing arguments to container: ${CONTAINER_ARGS[@]}"

# The container's output directory is the *value* of the --outputdir argument
# This maps the specified host directory to a directory with the same name inside the container.
docker run --rm \
    -v "${HOST_OUTPUT_PATH}:/app/${HOST_OUTPUT_DIR}" \
    ${IMAGE_NAME}:${TAG} conda run -n omics_cancer_explore python omics_cancer_explore/omics_cancer_explore.py "${CONTAINER_ARGS[@]}"

echo "Analysis complete. Check the '${HOST_OUTPUT_DIR}' directory for results."
