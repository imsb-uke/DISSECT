config = {
    "experiment_folder": "/home/user/expriment",  # Path to save outputs. Default: expriment

    "simulation_params": { 
        "scdata": "/home/user/experiment/data.h5ad",  # Path to sc/snRNA-seq data, should be anndata
        "n_samples": None,  # Number of samples to generate. Default: 1000 times the number of celltypes,
        "type": "bulk",
        "celltype_col": "celltype",  # Name of the column corresponding to cell-type labels in adata.obs
        "batch_col": None,  # If more than one batches are present, name of the column corrsponding to batch labels in adata.obs
        "cells_per_sample": None,  # Number of cells to sample to generate one sample.
        # Default 500
        "downsample": None,  # If simulation_type is ST, a float is used to downsample counts
        "preprocess": None,
        "filter": {  # Filtering of sc/snRNA-seq before simulating
            "min_genes": 200,
            "min_cells": 3,
            "mt_cutoff": 5,
            "min_expr": 0,  # in log2(1+count)
        },
        "concentration": None,  # Concentration parameter for dirichlet distribution
        # Should be a vector of same length as the number of cell-types with non-zero values
        # Higher concentrations will be favored. e.g. concentration [0.2,0.2,1] for 3 cell-types will make fractions
        # of the third cell-types higher.
        # Default: Vector of ones.
        "prop_sparse": 0.5,  # Proportion of sparse samples to generate. Default: 0.5
        # Sparse samples are samples in which some cell-types do not exist.
        # Probabilities of cell-types to not be present in the generate sample are uniform.
        "generate_component_figures": True,  # Computes PCA of celltype signatures per generated sample
    },

    "deconv_params": {
        "test_dataset": "../bulk.txt",
        "test_dataset_format": "txt",  # Either tab-delimited txt file with genes in rows or h5ad file compatible with Scanpy.
        "test_dataset_type": "bulk",  # bulk, microarray or spatial
        "duplicated": "first",  # In case, there are duplicated genes in the test_dataset. To use the first occuring gene, write first. To sum the duplicated genes, write sum. To take average, write mean
        "normalize_simulated": "cpm",  # "cpm", # Only CPM and None is supported. Write CPM if not already TPM/CPM.
        "normalize_test": "cpm",  # Write CPM if not already TPM/CPM
        "var_cutoff": 0.1,  # variance cutoff for gene filtering
        "test_in_mix": None,  # Number of test samples to use in the generation of online mixtures. None uses all samples.
        "simulated": True,  # True if dataset is already simulated. False, if it is a single-cell dataset.
        "sig_matrix": False,
        "mix": "srm",
        "save_config": True,
        "network_params": {
            "n_hidden_layers": 4,  # Number of hidden layers
            "hidden_units": [
                512,
                256,
                128,
                64,
                ],  # Sizes of the hidden dense layers. The length of this list should be same as n_hidden_layers above.
            "hidden_activation": "relu6",  # Activation of hidden layers. Choose ones supported in keras or relu6.
            "output_activation": "softmax",  # Activation of output layer.
            "loss": "kldivergence",  # Options - kldivergence, l2, l1. KL divergence will only work properly if output activation is softmax.
            "n_steps": 5000,  # Number of training steps
            "lr": 1e-5,  # Learning rate
            "batch_size": 64,  # best - 64 # batch size
            "dropout": None,  # If you would like dropoouts in the model, write a list with same number of elements as n_hidden_layers above corresponding to each dropout layer.
            # An example is [0.2,0.2,0,0.3,0.1,0.2]
            "n_steps_expr": 5000
            },  # Parameters to use to build network.
        "alpha_range": [
            0.1,
            0.9,
            ],  # Alpha parameter to create mixtures per batch is uniformly sampled between these two values
        "normalization_per_batch": "log1p-MinMax",  # normalization of batches. Only per log1p-MinMax or None are supported
        "models": [1, 2, 3, 4, 5],
        
    },
}