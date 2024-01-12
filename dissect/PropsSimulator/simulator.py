import os
import sys
import random
import numpy as np
import pandas as pd
import anndata as ad
from anndata import AnnData
import scanpy as sc
import shutil
from tqdm import tqdm
import json
import matplotlib.pyplot as plt


class Simulate(object):
    def __init__(self):
        pass
    
    def initialize(self, config):
        self.config = config
        self.sc_adata = sc.read(config["simulation_params"]["scdata"])
        if "sparse" in str(type(self.sc_adata.X)):
            self.sc_adata.X = np.array(self.sc_adata.X.todense())
        self.celltypes = np.sort(
            np.array(self.sc_adata.obs[config["simulation_params"]["celltype_col"]].unique())
        )
        self.n_celltypes = len(self.celltypes)

    def generate_props(self):
        if not self.config["simulation_params"]["n_samples"]:
            self.config["simulation_params"]["n_samples"] = 1000 * self.n_celltypes
        if not self.config["simulation_params"]["cells_per_sample"]:
            self.config["simulation_params"]["cells_per_sample"] = 1000
        self.n_sparse = int(self.config["simulation_params"]["n_samples"] * self.config["simulation_params"]["prop_sparse"])
        self.n_complete = self.config["simulation_params"]["n_samples"] - self.n_sparse

        ##### Complete
        if not self.config["simulation_params"]["concentration"]:
            self.config["simulation_params"]["concentration"] = np.ones(self.n_celltypes)
        self.props_complete = np.random.dirichlet(
            self.config["simulation_params"]["concentration"], self.n_complete
        ).astype(np.float32)
        self.min_prc = 1 / self.config["simulation_params"]["cells_per_sample"]
        self.props_complete[self.props_complete < self.min_prc] = 0
        # Re-normalize
        self.props_complete = (
            self.props_complete.T / self.props_complete.sum(axis=1)
        ).T
        self.cells_complete = self.props_complete * self.config["simulation_params"]["cells_per_sample"]
        self.cells_complete = np.round(self.cells_complete, 0).astype(int)
        # Update props to maintain summing to 1
        self.props_complete = (
            self.cells_complete.T / self.cells_complete.sum(axis=1)
        ).T.astype(np.float32)

        ##### Sparse
        keep_sparse = np.ones((self.n_sparse, self.n_celltypes), dtype=np.float32)

        no_keep = np.random.randint(1, self.n_celltypes, size=self.n_sparse)
        no_keeps = []
        for i in no_keep:
            no_keeps.append(
                np.random.choice(list(range(self.n_celltypes)), size=i, replace=False)
            )

        for i in range(len(no_keep)):
            keep_sparse[i, no_keeps[i]] = 1e-6

        self.props_sparse = np.zeros(
            (self.n_sparse, self.n_celltypes), dtype=np.float32
        )
        for i in range(self.n_sparse):
            self.props_sparse[i, :] = np.random.dirichlet(keep_sparse[i], 1)[0]

        self.props_sparse[self.props_sparse < self.min_prc] = 0
        self.props_sparse = (self.props_sparse.T / self.props_sparse.sum(axis=1)).T

        self.cells_sparse = self.props_sparse * self.config["simulation_params"]["cells_per_sample"]
        self.cells_sparse = np.round(self.cells_sparse, 0).astype(int)
        self.props_sparse = (
            self.cells_sparse.T / self.cells_sparse.sum(axis=1)
        ).T.astype(np.float32)


        if not os.path.exists(self.config["experiment_folder"]):
            os.mkdir(self.config["experiment_folder"])
        else:
            from datetime import datetime

            now = datetime.now()
            date_time = now.strftime("%m.%d.%Y_%H.%M.%S")
            name = "experiment_" + date_time
            print(
                "Specified experiment_folder {} already exists. Creating a new folder with name {}.".format(
                    self.config["experiment_folder"], name
                )
            )
            self.config["experiment_folder"] = name
            os.mkdir(self.config["experiment_folder"])

        self.simulation_folder = os.path.join(self.config["experiment_folder"], "simulation")
        if not os.path.exists(self.simulation_folder):
            os.mkdir(self.simulation_folder)
        else:
            sys.exit(f"folder {self.simulation_folder} already exists.")

        fig = plt.figure()
        ax = plt.boxplot(self.props_complete, labels=self.celltypes)  #
        plt.ylabel("Proportion")
        plt.title("Proportions of cell-types in generated samples")
        plt.savefig(
            os.path.join(self.simulation_folder, "boxplot_props_complete.pdf")
        )
        fig = plt.figure()
        ax = plt.boxplot(self.cells_complete, labels=self.celltypes)
        plt.ylabel("Count")
        plt.title("Counts of cell-types in generated samples")
        plt.savefig(
            os.path.join(self.simulation_folder, "boxplot_ncells_complete.pdf")
        )

        fig = plt.figure()
        ax = plt.boxplot(self.props_sparse, labels=self.celltypes)  #
        plt.ylabel("Proportion")
        plt.title("Proportions of cell-types in generated samples")
        plt.savefig(
            os.path.join(self.simulation_folder, "boxplot_props_sparse.pdf")
        )
        fig = plt.figure()
        ax = plt.boxplot(self.cells_sparse, labels=self.celltypes)
        plt.ylabel("Count")
        plt.title("Counts of cell-types in generated samples")
        plt.savefig(
            os.path.join(self.simulation_folder, "boxplot_ncells_sparse.pdf")
        )

        self.props = np.concatenate([self.props_complete, self.props_sparse], axis=0)
        self.cells = np.concatenate([self.cells_complete, self.cells_sparse], axis=0)

    def preprocess(self):
        self.sc_adata.var_names_make_unique()
        sc.pp.filter_cells(self.sc_adata, min_genes=self.config["simulation_params"]["filter"]["min_genes"])
        sc.pp.filter_genes(self.sc_adata, min_cells=self.config["simulation_params"]["filter"]["min_cells"])

        self.sc_adata.var["mt"] = self.sc_adata.var_names.str.startswith(
            "MT-"
        )  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(
            self.sc_adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )
        self.sc_adata = self.sc_adata[
            self.sc_adata.obs.pct_counts_mt < self.config["simulation_params"]["filter"]["mt_cutoff"], :
        ]
        sc.pp.normalize_per_cell(self.sc_adata)

        tmp = self.sc_adata.copy()
        sc.pp.log1p(tmp)
        sc.pp.highly_variable_genes(tmp)
        heg = tmp.var[tmp.var.means > self.config["simulation_params"]["filter"]["min_expr"]].index
        self.sc_adata = self.sc_adata[:, heg]
        del tmp
        del heg

    def simulate(self, save=True):
        adatas = {}
        genes = self.sc_adata.var_names
        for celltype in self.celltypes:
            adatas[celltype] = np.array(
                self.sc_adata[
                    self.sc_adata.obs[self.config["simulation_params"]["celltype_col"]] == celltype
                ].X
            )
        self.sc_adata = None

        adata = AnnData(
            np.zeros((self.config["simulation_params"]["n_samples"], len(genes))),
            dtype=np.float32,
            var=pd.DataFrame(index=genes),
            obs=pd.DataFrame(self.props, columns=self.celltypes),
        )
        adata.obsm["cells"] = self.cells
        for celltype in self.celltypes:
            adata.layers[celltype] = np.zeros(adata.X.shape, dtype=np.float32)
        for i in tqdm(range(self.cells.shape[0])):
            for j in range(self.n_celltypes):
                cell_idxs = np.random.choice(
                    range(adatas[self.celltypes[j]].shape[0]),
                    self.cells[i][j],
                    replace=True,
                )
                sample = adatas[self.celltypes[j]][cell_idxs, :].sum(axis=0)
                adata.layers[self.celltypes[j]][i, :] = sample
                if j == 0:
                    total_sample = sample
                else:
                    total_sample = total_sample + sample
            adata.X[i, :] = total_sample

        if save:
            adata.write(os.path.join(self.simulation_folder, "simulated.h5ad"))
            sc.set_figure_params(dpi=200)
            tmp = adata.copy()
            sc.pp.normalize_total(tmp, target_sum=1e6)
            sc.pp.log1p(tmp)
            sc.tl.pca(tmp)
            sc.pl.pca(tmp, color=adata.obs.columns, show=False)
            plt.savefig(
                os.path.join(
                    self.simulation_folder, "scatterplot_pca_simulated.pdf"
                )
            )

            if self.config["simulation_params"]["generate_component_figures"]:
                idxs = {}
                celltypes_col = []
                for j in range(self.n_celltypes):
                    celltype = self.celltypes[j]
                    idxs[celltype] = adata.obsm["cells"][:, j] > 0
                    celltypes_col = celltypes_col + [celltype] * idxs[celltype].sum()
                X = [
                    adata.layers[celltype][idxs[celltype], :]
                    for celltype in self.celltypes
                ]
                tmp1 = AnnData(
                    np.concatenate(X, axis=0),
                    obs=pd.DataFrame(celltypes_col, columns=["Celltype"]),
                )

                sc.pp.subsample(tmp1, fraction=0.5, random_state=42)
                sc.pp.normalize_total(tmp1, target_sum=1e6)
                sc.pp.log1p(tmp1)
                sc.tl.pca(tmp1)
                fig, ax = plt.subplots(nrows=1, ncols=1)
                sc.pl.pca(
                    tmp1, color="Celltype", show=False, legend_loc="on data", ax=ax
                )
                plt.savefig(
                    os.path.join(
                        self.simulation_folder,
                        "scatterplot_pca_simulated_celltypes.pdf",
                    ),
                    bbox_inches="tight",
                )
        else:
            return adata

    def simulate_per_batch(self, save=True):
        adata_orig = self.sc_adata.copy()
        adatas = []
        for batch in self.batches:
            self.sc_adata = adata_orig[
                adata_orig.obs[self.config["simulation_params"]["batch_col"]] == batch
            ]
            adatas.append(self.simulate(save=False))
        adata = adatas[0].concatenate(adatas[1:], batch_categories=self.batches)

        adata.write(os.path.join(self.simulation_folder, "simulated.h5ad"))

        if save:
            sc.set_figure_params(dpi=200)
            tmp = adata.copy()
            sc.pp.normalize_total(tmp, target_sum=1e6)
            sc.pp.log1p(tmp)
            sc.tl.pca(tmp)
            sc.pl.pca(tmp, color=adata.obs.columns, show=False)
            plt.savefig(
                os.path.join(
                    self.simulation_folder, "scatterplot_pca_simulated.pdf"
                )
            )

            adata.write(os.path.join(self.simulation_folder, "simulated.h5ad"))
            sc.set_figure_params(dpi=200)
            tmp = adata.copy()
            sc.pp.normalize_total(tmp, target_sum=1e6)
            sc.pp.log1p(tmp)
            sc.tl.pca(tmp)
            sc.pl.pca(tmp, color=adata.obs.columns, show=False)
            plt.savefig(
                os.path.join(
                    self.simulation_folder, "scatterplot_pca_simulated.pdf"
                )
            )

            if self.config["simulation_params"]["generate_component_figures"]:
                idxs = {}
                celltypes_col = []
                for j in range(self.n_celltypes):
                    celltype = self.celltypes[j]
                    idxs[celltype] = adata.obsm["cells"][:, j] > 0
                    celltypes_col = celltypes_col + [celltype] * idxs[celltype].sum()
                X = [
                    adata.layers[celltype][idxs[celltype], :]
                    for celltype in self.celltypes
                ]
                tmp1 = AnnData(
                    np.concatenate(X, axis=0),
                    obs=pd.DataFrame(celltypes_col, columns=["Celltype"]),
                )

                sc.pp.subsample(tmp1, fraction=0.5, random_state=42)
                sc.pp.normalize_total(tmp1, target_sum=1e6)
                sc.pp.log1p(tmp1)
                sc.tl.pca(tmp1)
                fig, ax = plt.subplots(nrows=1, ncols=1)
                sc.pl.pca(
                    tmp1, color="Celltype", show=False, legend_loc="on data", ax=ax
                )
                plt.savefig(
                    os.path.join(
                        self.simulation_folder,
                        "scatterplot_pca_simulated_celltypes.pdf",
                    ),
                    bbox_inches="tight",
                )

class Simulate_st(object):
    def __init__(self):
        pass

    def initialize(self, config):
        self.config["simulation_params"] = config
        self.sc_adata = sc.read(config["simulation_params"]["scdata"])
        if "sparse" in str(type(self.sc_adata.X)):
            self.sc_adata.X = np.array(self.sc_adata.X.todense())
        self.celltypes = np.sort(
            np.array(self.sc_adata.obs[config["simulation_params"]["celltype_col"]].unique())
        )
        self.n_celltypes = len(self.celltypes)

    def generate_props(self):
        if not self.config["simulation_params"]["n_samples"]:
            self.config["simulation_params"]["n_samples"] = 1000 * self.n_celltypes
        self.n_sparse = self.config["simulation_params"]["n_samples"]
        ##### Sparse
        keep_sparse = np.ones((self.n_sparse, self.n_celltypes), dtype=np.float32)

        keep = np.random.randint(1, 6, size=self.n_sparse)
        keeps = []
        for i in keep:
            keeps.append(
                np.random.choice(list(range(self.n_celltypes)), size=i, replace=False)
            )

        for i in range(len(keep)):
            keep_sparse[i, ~keeps[i]] = 1e-6

        self.props_sparse = np.zeros(
            (self.n_sparse, self.n_celltypes), dtype=np.float32
        )
        for i in range(self.n_sparse):
            self.props_sparse[i, :] = np.random.dirichlet(keep_sparse[i], 1)[0]

        self.props_sparse = (self.props_sparse.T / self.props_sparse.sum(axis=1)).T

        for i in range(self.props_sparse.shape[0]):
            n = np.random.choice(list(range(5, 12)), size=1, replace=False)
            self.cells_sparse = self.props_sparse * n
        self.cells_sparse = np.round(self.cells_sparse, 0).astype(int)
        self.props_sparse = (
            self.cells_sparse.T / self.cells_sparse.sum(axis=1)
        ).T.astype(np.float32)

        if not os.path.exists(self.config["experiment_folder"]):
            os.mkdir(self.config["experiment_folder"])
        else:
            from datetime import datetime

            now = datetime.now()
            date_time = now.strftime("%m.%d.%Y_%H.%M.%S")
            name = "experiment_" + date_time
            print(
                "Specified experiment_folder {} already exists. Creating a new folder with name {}.".format(
                    self.config["experiment_folder"], name
                )
            )
            self.config["experiment_folder"] = name
            os.mkdir(self.config["experiment_folder"])

        self.simulation_folder = os.path.join(self.config["experiment_folder"], "simulation")
        if not os.path.exists(self.simulation_folder):
            os.mkdir(self.simulation_folder)
        else:
            sys.exit(f"folder {self.simulation_folder} already exists.")

        fig = plt.figure()
        ax = plt.boxplot(self.props_sparse, labels=self.celltypes)  #
        plt.ylabel("Proportion")
        plt.title("Proportions of cell-types in generated samples")
        plt.savefig(
            os.path.join(self.config.simulation_folder, "boxplot_props_sparse.pdf")
        )
        fig = plt.figure()
        ax = plt.boxplot(self.cells_sparse, labels=self.celltypes)
        plt.ylabel("Count")
        plt.title("Counts of cell-types in generated samples")
        plt.savefig(
            os.path.join(self.config.simulation_folder, "boxplot_ncells_sparse.pdf")
        )

        self.props = self.props_sparse
        self.cells = self.cells_sparse

    def preprocess(self):
        self.sc_adata.var_names_make_unique()
        sc.pp.filter_cells(self.sc_adata, min_genes=self.config["simulation_params"]["filter"]["min_genes"])
        sc.pp.filter_genes(self.sc_adata, min_cells=self.config["simulation_params"]["filter"]["min_cells"])

        self.sc_adata.var["mt"] = self.sc_adata.var_names.str.startswith(
            "MT-"
        )  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(
            self.sc_adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )
        self.sc_adata = self.sc_adata[
            self.sc_adata.obs.pct_counts_mt < self.config["simulation_params"]["filter"]["mt_cutoff"], :
        ]
        sc.pp.normalize_per_cell(self.sc_adata)

        tmp = self.sc_adata.copy()
        sc.pp.log1p(tmp)
        sc.pp.highly_variable_genes(tmp)
        heg = tmp.var[tmp.var.means > self.config["simulation_params"]["filter"]["min_expr"]].index
        self.sc_adata = self.sc_adata[:, heg]
        del tmp
        del heg

    def simulate(self, save=True):
        adatas = {}
        genes = self.sc_adata.var_names
        for celltype in self.celltypes:
            adatas[celltype] = np.array(
                self.sc_adata[
                    self.sc_adata.obs[self.config["simulation_params"]["celltype_col"]] == celltype
                ].X
            )
        self.sc_adata = None

        adata = AnnData(
            np.zeros((self.config["simulation_params"]["n_samples"], len(genes))),
            var=pd.DataFrame(index=genes),
            obs=pd.DataFrame(self.props, columns=self.celltypes),
        )
        adata.obsm["cells"] = self.cells
        for celltype in self.celltypes:
            adata.layers[celltype] = np.zeros(adata.X.shape, dtype=np.float32)
        for i in tqdm(range(self.cells.shape[0])):
            for j in range(self.n_celltypes):
                cell_idxs = np.random.choice(
                    range(adatas[self.celltypes[j]].shape[0]),
                    self.cells[i][j],
                    replace=True,
                )
                sample = adatas[self.celltypes[j]][cell_idxs, :].sum(axis=0)
                adata.layers[self.celltypes[j]][i, :] = sample
                if j == 0:
                    total_sample = sample
                else:
                    total_sample = total_sample + sample
            adata.X[i, :] = total_sample

        if self.config["simulation_params"]["downsample"]:
            sc.pp.downsample_counts(
                adata, counts_per_cell=self.config["simulation_params"]["downsample"] * adata.X.sum(1)
            )
        if save:
            adata.write(os.path.join(self.simulation_folder, "simulated.h5ad"))
            sc.set_figure_params(dpi=200)
            tmp = adata.copy()
            sc.pp.normalize_total(tmp, target_sum=1e6)
            sc.pp.log1p(tmp)
            sc.tl.pca(tmp)
            sc.pl.pca(tmp, color=adata.obs.columns, show=False)
            plt.savefig(
                os.path.join(
                    self.simulation_folder, "scatterplot_pca_simulated.pdf"
                )
            )

            if self.config["generate_component_figures"]:
                idxs = {}
                celltypes_col = []
                for j in range(self.n_celltypes):
                    celltype = self.celltypes[j]
                    idxs[celltype] = adata.obsm["cells"][:, j] > 0
                    celltypes_col = celltypes_col + [celltype] * idxs[celltype].sum()
                X = [
                    adata.layers[celltype][idxs[celltype], :]
                    for celltype in self.celltypes
                ]
                tmp1 = AnnData(
                    np.concatenate(X, axis=0),
                    obs=pd.DataFrame(celltypes_col, columns=["Celltype"]),
                )

                sc.pp.subsample(tmp1, fraction=0.5, random_state=42)
                sc.pp.normalize_total(tmp1, target_sum=1e6)
                sc.pp.log1p(tmp1)
                sc.tl.pca(tmp1)
                fig, ax = plt.subplots(nrows=1, ncols=1)
                sc.pl.pca(
                    tmp1, color="Celltype", show=False, legend_loc="on data", ax=ax
                )
                plt.savefig(
                    os.path.join(
                        self.simulation_folder,
                        "scatterplot_pca_simulated_celltypes.pdf",
                    ),
                    bbox_inches="tight",
                )
        else:
            return adata

    def simulate_per_batch(self, save=True):
        adata_orig = self.sc_adata.copy()
        adatas = []
        for batch in self.batches:
            self.sc_adata = adata_orig[
                adata_orig.obs[self.config["simulation_params"]["batch_col"]] == batch
            ]
            adatas.append(self.simulate(save=False))
        adata = adatas[0].concatenate(adatas[1:], batch_categories=self.batches)

        adata.write(os.path.join(self.simulation_folder, "simulated.h5ad"))

        if save:
            sc.set_figure_params(dpi=200)
            tmp = adata.copy()
            sc.pp.normalize_total(tmp, target_sum=1e6)
            sc.pp.log1p(tmp)
            sc.tl.pca(tmp)
            sc.pl.pca(tmp, color=adata.obs.columns, show=False)
            plt.savefig(
                os.path.join(
                    self.simulation_folder, "scatterplot_pca_simulated.pdf"
                )
            )

            adata.write(os.path.join(self.simulation_folder, "simulated.h5ad"))
            sc.set_figure_params(dpi=200)
            tmp = adata.copy()
            sc.pp.normalize_total(tmp, target_sum=1e6)
            sc.pp.log1p(tmp)
            sc.tl.pca(tmp)
            sc.pl.pca(tmp, color=adata.obs.columns, show=False)
            plt.savefig(
                os.path.join(
                    self.simulation_folder, "scatterplot_pca_simulated.pdf"
                )
            )

            if self.config["simulation_params"]["generate_component_figures"]:
                idxs = {}
                celltypes_col = []
                for j in range(self.n_celltypes):
                    celltype = self.celltypes[j]
                    idxs[celltype] = adata.obsm["cells"][:, j] > 0
                    celltypes_col = celltypes_col + [celltype] * idxs[celltype].sum()
                X = [
                    adata.layers[celltype][idxs[celltype], :]
                    for celltype in self.celltypes
                ]
                tmp1 = AnnData(
                    np.concatenate(X, axis=0),
                    obs=pd.DataFrame(celltypes_col, columns=["Celltype"]),
                )

                sc.pp.subsample(tmp1, fraction=0.5, random_state=42)
                sc.pp.normalize_total(tmp1, target_sum=1e6)
                sc.pp.log1p(tmp1)
                sc.tl.pca(tmp1)
                fig, ax = plt.subplots(nrows=1, ncols=1)
                sc.pl.pca(
                    tmp1, color="Celltype", show=False, legend_loc="on data", ax=ax
                )
                plt.savefig(
                    os.path.join(
                        self.simulation_folder,
                        "scatterplot_pca_simulated_celltypes.pdf",
                    ),
                    bbox_inches="tight",
                )

def save_dict_to_file(config):
    f = open(os.path.join(config["simulation_params"]["simulation_folder"], "simulation_config.py"), "w")
    f.write("config = ")
    f.write(str(config))
    f.close()

def simulate(config):

    s = 42
    random.seed(s)
    np.random.seed(s)
    if config["simulation_params"]["type"]=="bulk":
        sim = Simulate()
    else:
        sim = Simulate_st()
    sim.initialize(config)
    sim.preprocess()
    sim.generate_props()
    batch_col = sim.config["simulation_params"]["batch_col"]
    columns = sim.sc_adata.obs.columns
    print(batch_col)
    print(columns)
    if batch_col in columns:
        sim.batches = np.array(sim.sc_adata.obs[sim.config["simulation_params"]["batch_col"]].unique())
        if len(sim.batches) > 1:
            sim.simulate_per_batch(save=True)
        else:
            sim.simulate(save=True)
    else:
        print("Number of batches in single-cell data is 1.")
        sim.simulate(save=True)
    sim.config["simulation_params"]["simulation_folder"] = os.path.join(sim.config["experiment_folder"], "simulation")
    sim.config["simulation_params"]["concentration"] = list(sim.config["simulation_params"]["concentration"])
    sim.config["deconv_params"]["reference"] = os.path.join(sim.config["simulation_params"]["simulation_folder"], "simulated.h5ad")
    save_dict_to_file(sim.config)
    