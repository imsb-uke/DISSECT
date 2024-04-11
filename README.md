# DISSECT

DISSECT is a deep semi-supervised learning framework that performs complete deconvolution of gene expression mixtures such as bulk RNAseq, proteomics and spatial transcriptomics

This repository contains the functions to run DISSECT interactively in Python. The code to reproduce figures and run deconvolution pipeline through comamnd line is available [here](https://github.com/imsb-uke/deconvolution).  

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. Folder tutorials include notebooks to get started.

## Prerequisites

conda >= v22 through [Anaconda](https://docs.anaconda.com/free/anaconda/install/index.html) or [miniconda](https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html).

## Installing
```shell
## Installation

# Create and activate virtual environment. This is recommended to avoid conflict in dependencies.
conda create -y -n dissect python=3.9
conda activate dissect

# Clone DISSECT
git clone https://github.com/imsb-uke/DISSECT

# Install dependencies and DISSECT

pip install -r DISSECT/requirements.txt --use-deprecated=legacy-resolver
pip install DISSECT/.

# Install jupyter lab for interactive execution of notebooks
conda install -y jupyter

## Tutorials to get started
# Go to tutorials directory within DISSECT
cd DISSECT/tutorials

## Launch jupyter lab
jupyter notebook

```
## GPU usage
By default, tensorflow-gpu which is installed while installed DISSECT works as long as appropriate CUDA driver is installed. DISSECT uses tensorflow-gpu version 2.7.0 with CUDA 11.2 and cuDNN 8.1. The available devices to tensorflow can be checked as below.

```python
import tensorflow as tf
gpus = tf.config.list_physical_devices("GPU")
print(gpus)

```
This will output a list of the available GPU devices as the output below where we have 1 GPU available. 
```
[PhysicalDevice(name='/physical_device:GPU:0', device_type='GPU')] 
```
In case there are multiple GPUs available, a particular GPU can be set by,

```python
gpu_number = 0 # Using only the first GPU
tf.config.experimental.set_visible_devices(gpus[gpu_number], 'GPU')
```

## Tutorials
Interactive tutorials including required data are available as part of this repository at [Tutorials](https://github.com/imsb-uke/DISSECT/tree/main/tutorials).
1. [tutorial.ipynb](https://github.com/imsb-uke/DISSECT/tree/main/tutorials/tutorial.ipynb): Step by step deconvolution for bulk
2. [tutorial.ipynb](https://github.com/imsb-uke/DISSECT/tree/main/tutorials/tutorial_spatial.ipynb): Step by step deconvolution of spatial transcriptomics data (10x Visium)

### A note on the supported formats
bulk RNAseq/Proteomics: A tab seperated file (gene symbols in rows, sample names in columns), or an h5ad compatible with Scanpy >=1.8.0.
scRNAseq: An h5ad compatible with Scanpy >=1.8.0. The cell metadata (.obs) must contain cell type labels (By default in column cell_type)
Spatial transcriptomics: A h5ad compatible with Scanpy >=1.8.0.
Gene ID formats should match between scRNAseq and bulk/spatial data but can be anything (HGNC symbols, ENSEMBL etc.).



To get answers quickly for a problem or feature request, please open an [issue](https://github.com/imsb-uke/DISSECT/issues).
