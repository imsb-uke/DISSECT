# DISSECT

DISSECT is a deep semi-supervised learning framework that performs complete deconvolution of gene expression mixtures such as bulk RNAseq, proteomics and spatial transcriptomics

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. Folder tutorials include notebooks to get started.

## Prerequisites

conda >= v22 through [Anaconda](https://docs.anaconda.com/free/anaconda/install/index.html) or [miniconda](https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html).


## Installing
```python
## Installation

# Create and activate virtual environment
conda create -y -n dissect python=3.8
conda activate dissect

# Clone DISSECT
git clone https://github.com/imsb-uke/DISSECT

# Install dependencies and DISSECT

pip install -r DISSECT/requirements.txt --use-deprecated=legacy-resolver
pip install DISSECT/.

# Install jupyter lab
conda install -y jupyter

## Tutorials to get started
# Go to tutorials directory within DISSECT
cd DISSECT/tutorials

```

## Tutorials
Interactive tutorials including required data are available as part of this repository at [Tutorials](https://github.com/imsb-uke/DISSECT/tree/main/tutorials).
1. expanded_tutorial.ipynb: Step by step deconvolution
2. minimal_tutorial.ipynb: Complete deconvolution from a single configuration file using minimal steps of code

To get answers quickly for a problem or feature request, please open an issue.
