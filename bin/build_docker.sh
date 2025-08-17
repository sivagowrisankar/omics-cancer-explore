#!/bin/bash
set -e

IMAGE_NAME="omics_cancer_analyzer"
TAG="latest"

echo "Building Docker image: ${IMAGE_NAME}:${TAG}"

docker build -t ${IMAGE_NAME}:${TAG} -f docker/omics_cancer_explore/Dockerfile .

echo "Docker image built successfully."
