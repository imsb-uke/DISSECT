import os
import sys
import numpy as np
import pandas as pd
from anndata import AnnData
import scanpy as sc
import shutil
import json

def dataset(config):
    """
    Prepares datasets for training dissect. Uses parameters from config.py.
    """
    savedir = config["experiment_folder"]
    dataset_path = os.path.join(savedir, "datasets")
    simulation_path = os.path.join(savedir, "simulation")

    if not os.path.exists(savedir):
        sys.exit("Directory doesn't exist. Please run simulation before.")
    elif os.path.exists(savedir) and not os.path.exists(simulation_path):
        sys.exit("Simulations do not exist. Please run simulation before.")
    elif os.path.exists(savedir) and os.path.exists(simulation_path):
        if os.path.exists(dataset_path):
            sys.exit("Datasets are already prepared for the given configuration.")
        else:
            os.mkdir(dataset_path)

    if config["deconv_params"]["save_config"]:
        with open(
            os.path.join(config["experiment_folder"], "main_config.json"), "w"
        ) as f:
            json.dump(config, f)
    # if not os.path.exists(os.path.join(savedir, "main_config.py")):
    #    shutil.copyfile("configs/main_config.py", os.path.join(savedir, "main_config.py"))
    ###################
    # Read test dataset
    ###################
    if config["deconv_params"]["test_dataset_format"] == "txt":
        X_real = pd.read_table(config["deconv_params"]["test_dataset"], index_col=0)
    elif config["deconv_params"]["test_dataset_format"] == "h5ad":
        X_real = sc.read(config["deconv_params"]["test_dataset"])
        X_real = pd.DataFrame(
            X_real.X, index=X_real.obs.index.tolist(), columns=X_real.var_names.tolist()
        ).T

    if config["deconv_params"]["test_dataset_type"] == "microarray":
        X_real = 2**X_real - 1
        X_real[X_real < 0] = 0

    # Filter genes
    if config["deconv_params"]["var_cutoff"]:
        print(
            "Removing genes which have less than {} variance in their expressions.".format(
                config["deconv_params"]["var_cutoff"]
            )
        )
        genes_to_keep = X_real.var(axis=1) > config["deconv_params"]["var_cutoff"]
        X_real = X_real.loc[genes_to_keep]

    X_real = X_real.T  # s,g
    genes_real = X_real.columns.tolist()
    sample_names = X_real.index.tolist()

    n_genes = len(genes_real)
    n_distinct_genes = len(set(genes_real))

    if n_genes != n_distinct_genes:
        if config["deconv_params"]["duplicated"] == "first":
            X_real = X_real.loc[:, ~X_real.columns.duplicated(keep="first")]
            s = "Kept the gene expressions of the first occured gene"
            print(
                "There are duplicated genes in the test dataset. {} as specified in the parameter duplicates in main_config.py".format(
                    s
                )
            )
        elif config["deconv_params"]["duplicated"] == "sum":
            X_real = X_real.groupby(X_real.columns, axis=1).sum()
            s = "Summed the genes expressions of duplicated gene names"
            print(
                "There are duplicated genes in the test dataset. {} as specified in the parameter duplicates in main_config.py".format(
                    s
                )
            )
        elif config["deconv_params"]["duplicated"] == "sum":
            X_real = X_real.groupby(X_real.columns, axis=1).mean()
            s = "Took arithmetic mean of the genes expressions of duplicated gene names"
            print(
                "There are duplicated genes in the test dataset. {} as specified in the parameter duplicates in main_config.py".format(
                    s
                )
            )
        else:
            sys.exit(
                "duplicated setting {} not supported. Please check config.py file.".format(
                    config["deconv_params"]["duplicated"]
                )
            )

    print("test dataset has {} distinct and variable genes.".format(n_distinct_genes))

    # Normalization
    if config["deconv_params"]["normalize_test"] == "cpm":
        X_real_anndata = AnnData(
            np.array(X_real), dtype=np.float32, var=pd.DataFrame(index=X_real.columns.tolist())
        )
        sc.pp.normalize_total(X_real_anndata, target_sum=1e6)
        X_real = pd.DataFrame(
            X_real_anndata.X,
            index=X_real.index.tolist(),
            columns=X_real.columns.tolist(),
        )
        del X_real_anndata
    elif not config["deconv_params"]["normalize_test"]:
        pass
    else:
        sys.exit(
            "{} in normalize_test in config is not supported.".format(
                config["deconv_params"]["normalize_simulated"]
            )
        )

    ################
    # Read Reference
    ################
    X_sc = sc.read(config["deconv_params"]["reference"])

    # Simulated if not simulated
    if config["deconv_params"]["simulated"]:
        X_sim = X_sc
    else:
        X_sim = simulate(X_sc, config["deconv_params"]["simulation_params"])

    # Normalization
    if config["deconv_params"]["normalize_simulated"] == "cpm":
        sc.pp.normalize_total(X_sim, target_sum=1e6)
    elif not config["deconv_params"]["normalize_simulated"]:
        pass
    else:
        sys.exit(
            "{} in normalize_simulated in config is not supported.".format(
                config["deconv_params"]["normalize_simulated"]
            )
        )
    y_sim = X_sim.obs[
        [col for col in X_sim.obs.columns if col not in ["ds", "batch"]]
    ]  # X_sim.uns["cell_types"].tolist()
    y_sim[y_sim < 0.005] = 0
    y_sim = y_sim.div(y_sim.sum(1), 0)
    X_sim = pd.DataFrame(
        X_sim.X, columns=X_sim.var_names.tolist(), index=X_sim.obs.index.tolist()
    )

    genes_sim = X_sim.columns.tolist()
    n_genes = len(genes_sim)
    n_distinct_genes = len(set(genes_sim))

    if n_genes != n_distinct_genes:
        if config["deconv_params"]["duplicated"] == "first":
            X_sim = X_sim.loc[:, ~X_sim.columns.duplicated(keep="first")]
            s = "Kept the gene expressions of the first occured gene"
            print(
                "There are duplicated genes in the simulated dataset. {} as specified in the parameter duplicates in main_config.py".format(
                    s
                )
            )
        elif config["deconv_params"]["duplicated"] == "sum":
            X_sim = X_sim.groupby(X_sim.columns, axis=1).sum()
            s = "Summed the genes expressions of duplicated gene names"
            print(
                "There are duplicated genes in the simulated dataset. {} as specified in the parameter duplicates in main_config.py".format(
                    s
                )
            )
        elif config["deconv_params"]["duplicated"] == "sum":
            X_sim = X_sim.groupby(X_sim.columns, axis=1).mean()
            s = "Took arithmetic mean of the genes expressions of duplicated gene names"
            print(
                "There are duplicated genes in the simulated dataset. {} as specified in the parameter duplicates in main_config.py".format(
                    s
                )
            )
        else:
            sys.exit(
                "duplicated setting {} not supported. Please check config.py file.".format(
                    config["deconv_params"]["duplicated"]
                )
            )
    print("simulated dataset has {} distinct genes.".format(n_distinct_genes))

    # Prepare datasets
    genes_intersect = list(set(genes_real) & set(genes_sim))
    print(
        "There are {} common genes between simulated and test dataset.".format(
            len(genes_intersect)
        )
    )
    X_sim, X_real = X_sim.loc[:, genes_intersect], X_real.loc[:, genes_intersect]

    X_real_test = X_real.copy()
    if config["deconv_params"]["test_in_mix"]:
        X_real = X_real.iloc[0 : config["deconv_params"]["test_in_mix"], :]

    real_size = X_real.shape[0]
    sim_size = X_sim.shape[0]

    if X_real.shape[0] < X_sim.shape[0]:
        X_real = pd.concat([X_real] * int(sim_size / real_size), ignore_index=True)
    real_size = X_real.shape[0]
    curr_size = real_size

    while curr_size < sim_size:
        X_real = X_real.append(X_real.iloc[curr_size - real_size, :], ignore_index=True)
        curr_size += 1

    if X_sim.shape[0] < X_real.shape[0]:
        X_sim = pd.concat([X_sim] * int(real_size / sim_size), ignore_index=True)
    sim_size = X_sim.shape[0]
    curr_size = sim_size

    while curr_size < real_size:
        X_sim = X_sim.append(X_sim.iloc[curr_size - sim_size, :], ignore_index=True)
        curr_size += 1

    final_dataset = [X_real, X_sim, y_sim, X_real_test, sample_names]

    print("Saving numpy files.")

    np.save(os.path.join(dataset_path, "X_real_train.npy"), X_real, allow_pickle=True)
    np.save(os.path.join(dataset_path, "X_sim.npy"), X_sim, allow_pickle=True)
    np.save(os.path.join(dataset_path, "y_sim.npy"), y_sim, allow_pickle=True)
    np.save(
        os.path.join(dataset_path, "X_real_test.npy"), X_real_test, allow_pickle=True
    )
    celltypes = pd.DataFrame(index=y_sim.columns.tolist())
    celltypes.to_csv(os.path.join(dataset_path, "celltypes.txt"), sep="\t")
    sample_names = pd.DataFrame(index=sample_names)
    sample_names.to_csv(os.path.join(dataset_path, "sample_names.txt"), sep="\t")
    genes = pd.DataFrame(index=genes_intersect)
    genes.to_csv(os.path.join(dataset_path, "genes.txt"), sep="\t")

    print("Done.")


if __name__ == "__main__":
    from configs.main_config import config
    dataset(config)