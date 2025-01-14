import os
import sys
import random
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse
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
        self.save_expr = config["simulation_params"]["save_expr"]
        self.sc_adata.obs[config["simulation_params"]["celltype_col"]] = self.sc_adata.obs[config["simulation_params"]["celltype_col"]].astype(str)
        self.sc_adata.obs[config["simulation_params"]["celltype_col"]].replace("/","_",regex=True,inplace=True)
        self.sc_adata.obs[config["simulation_params"]["celltype_col"]] = self.sc_adata.obs[config["simulation_params"]["celltype_col"]].astype("category")
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
            self.config["simulation_params"]["cells_per_sample"] = 100
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
        if self.n_celltypes>10:
            plt.xticks(rotation=45, ha="right")
        plt.ylabel("Proportion")
        plt.title("Proportions of cell-types in generated samples")
        plt.savefig(
            os.path.join(self.simulation_folder, "boxplot_props_complete.pdf"),
            bbox_inches="tight"
        )
        # fig = plt.figure()
        # ax = plt.boxplot(self.cells_complete, labels=self.celltypes)
        # plt.ylabel("Count")
        # plt.title("Counts of cell-types in generated samples")
        # plt.savefig(
        #     os.path.join(self.simulation_folder, "boxplot_ncells_complete.pdf")
        # )

        # fig = plt.figure()
        # ax = plt.boxplot(self.props_sparse, labels=self.celltypes)  #
        # plt.ylabel("Proportion")
        # plt.title("Proportions of cell-types in generated samples")
        # plt.savefig(
        #     os.path.join(self.simulation_folder, "boxplot_props_sparse.pdf")
        # )
        # fig = plt.figure()
        # ax = plt.boxplot(self.cells_sparse, labels=self.celltypes)
        # plt.ylabel("Count")
        # plt.title("Counts of cell-types in generated samples")
        # plt.savefig(
        #     os.path.join(self.simulation_folder, "boxplot_ncells_sparse.pdf")
        # )

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
        if self.save_expr:
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
                if self.save_expr:
                    adata.layers[self.celltypes[j]][i, :] = sample
                if j == 0:
                    total_sample = sample
                else:
                    total_sample = total_sample + sample
            adata.X[i, :] = total_sample
        
        # convert to sparse matrices
        if self.save_expr:
            for i in range(self.n_celltypes):
                adata.layers[self.celltypes[i]] = csr_matrix(adata.layers[self.celltypes[i]])
        if save:
            adata.write(os.path.join(self.simulation_folder, "simulated.h5ad"))

            if self.config["simulation_params"]["generate_component_figures"]:
                sc.set_figure_params(dpi=200)
                tmp = adata.copy()
                sc.pp.normalize_total(tmp, target_sum=1e6)
                sc.pp.log1p(tmp)
                sc.tl.pca(tmp)
                sc.pl.pca(tmp, color=adata.obs.columns,show=False)
                plt.savefig(
                    os.path.join(
                        self.simulation_folder, "scatterplot_pca_simulated.pdf"
                    )
                )
                idxs = {}
                celltypes_col = []
                for j in range(self.n_celltypes):
                    celltype = self.celltypes[j]
                    idxs[celltype] = adata.obsm["cells"][:, j] > 0
                    celltypes_col = celltypes_col + [celltype] * idxs[celltype].sum()
                X = [
                    adata.layers[celltype].toarray()[idxs[celltype], :]
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
                    tmp1, color="Celltype", show=False, ax=ax
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

        # adata.write(os.path.join(self.simulation_folder, "simulated.h5ad"))

        if save:
            adata.write(os.path.join(self.simulation_folder, "simulated.h5ad"))
            
            if self.config["simulation_params"]["generate_component_figures"]:
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
                idxs = {}
                celltypes_col = []
                for j in range(self.n_celltypes):
                    celltype = self.celltypes[j]
                    idxs[celltype] = adata.obsm["cells"][:, j] > 0
                    celltypes_col = celltypes_col + [celltype] * idxs[celltype].sum()
                X = [
                    adata.layers[celltype].toarray()[idxs[celltype], :]
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
                    tmp1, color="Celltype", show=False, ax=ax
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
        self.config = config
        self.sc_adata = sc.read(config["simulation_params"]["scdata"])
        self.save_expr = config["simulation_params"]["save_expr"]
        self.sc_adata.obs[config["simulation_params"]["celltype_col"]] = self.sc_adata.obs[config["simulation_params"]["celltype_col"]].astype(str)
        self.sc_adata.obs[config["simulation_params"]["celltype_col"]].replace("/","_",regex=True,inplace=True)
        self.sc_adata.obs[config["simulation_params"]["celltype_col"]] = self.sc_adata.obs[config["simulation_params"]["celltype_col"]].astype("category")
        self.celltypes = np.sort(
            np.array(self.sc_adata.obs[config["simulation_params"]["celltype_col"]].unique())
        )
        self.n_celltypes = len(self.celltypes)
        
        self.celltype_indices = {}
        for ct in self.celltypes:
            mask = self.sc_adata.obs[config["simulation_params"]["celltype_col"]] == ct
            self.celltype_indices[ct] = np.where(mask)[0]

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
            idx_exclude = np.ones(self.n_celltypes, dtype=bool)
            idx_exclude[keeps[i]] = False
            keep_sparse[i, idx_exclude] = 1e-6

        # Generate proportions using Dirichlet
        props = np.zeros((self.n_sparse, self.n_celltypes), dtype=np.float32)
        for i in range(self.n_sparse):
            props[i, :] = np.random.dirichlet(keep_sparse[i], 1)[0]

        cells = np.zeros_like(props, dtype=int)
        for i in range(props.shape[0]):
            n_cells = np.random.choice(list(range(5, 12)), size=1, replace=False)[0]
            cells[i] = np.round(props[i] * n_cells).astype(int)

        # Recalculate proportions
        self.props = (cells.T / cells.sum(axis=1)).T.astype(np.float32)
        self.cells = cells

        self._setup_folders()
        self._create_plots()

    def _setup_folders(self):
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

    def _create_plots(self):
        fig = plt.figure()
        ax = plt.boxplot(self.props, labels=self.celltypes)
        if self.n_celltypes > 10:
            plt.xticks(rotation=45, ha="right")
        plt.ylabel("Proportion")
        plt.title("Proportions of cell-types in generated samples")
        plt.savefig(
            os.path.join(self.simulation_folder, "boxplot_props_sparse.pdf"),
            bbox_inches="tight"
        )

    def preprocess(self):
        self.sc_adata.var_names_make_unique()
        sc.pp.filter_cells(self.sc_adata, min_genes=self.config["simulation_params"]["filter"]["min_genes"])
        sc.pp.filter_genes(self.sc_adata, min_cells=self.config["simulation_params"]["filter"]["min_cells"])

        self.sc_adata.var["mt"] = self.sc_adata.var_names.str.startswith(
            "MT-"
        )
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
        genes = self.sc_adata.var_names
        n_cells_total = self.sc_adata.shape[0]  
        n_samples = self.cells.shape[0]         
        
        a = csr_matrix((n_cells_total, n_samples), dtype=np.float32)
        
        if self.save_expr:
            ct_arr = {}
            
        for sample_num in tqdm(range(n_samples)):
            if self.save_expr:
                ct_arr[sample_num] = {}
                
            row_indices = []
            col_indices = []
            data_values = []
            
            for ct_idx, ct in enumerate(self.celltypes):
                cell_num = self.cells[sample_num, ct_idx]
                if cell_num > 0:
                    valid_indices = self.celltype_indices[ct]
                    
                    cell_num = min(cell_num, len(valid_indices))
                    
                    if cell_num > 0:
                        sampled_indices = np.random.choice(valid_indices, cell_num, replace=True)
                        row_indices.extend(sampled_indices)
                        col_indices.extend([sample_num] * cell_num)
                        data_values.extend([1.0] * cell_num)
                        
                        if self.save_expr:
                            ct_arr[sample_num][ct] = self.sc_adata[sampled_indices].to_df().mean(0)
            
            if row_indices:
                max_idx = max(row_indices)
                if max_idx >= n_cells_total:
                    print(f"Warning: Found index {max_idx} exceeding {n_cells_total} at sample {sample_num}")
                    continue
                    
                new_data = csr_matrix(
                    (data_values, (row_indices, col_indices)),
                    shape=(n_cells_total, n_samples)
                )
                a = a + new_data

        S = self.sc_adata.X.T * a
        
        adata = AnnData(
            S.T,
            var=pd.DataFrame(index=genes),
            obs=pd.DataFrame(self.props, columns=self.celltypes),
        )
        adata.obsm["cells"] = self.cells
        
        if self.save_expr:
            for ct in self.celltypes:
                layer_data = np.zeros((adata.shape[0], adata.shape[1]), dtype=np.float32)
                for i in ct_arr:
                    if ct in ct_arr[i]:
                        layer_data[i] = ct_arr[i][ct]
                adata.layers[ct] = csr_matrix(layer_data)

        if self.config["simulation_params"]["downsample"]:
            sc.pp.downsample_counts(
                adata, counts_per_cell=self.config["simulation_params"]["downsample"] * adata.X.sum(1)
            )

        if save:
            adata.write(os.path.join(self.simulation_folder, "simulated.h5ad"))
            if self.config["simulation_params"]["generate_component_figures"]:
                self._generate_component_figures(adata)
        else:
            return adata

    def _generate_component_figures(self, adata):
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

        idxs = {}
        celltypes_col = []
        for j, celltype in enumerate(self.celltypes):
            idxs[celltype] = adata.obsm["cells"][:, j] > 0
            celltypes_col.extend([celltype] * idxs[celltype].sum())

        X = [adata.layers[celltype].toarray()[idxs[celltype], :] for celltype in self.celltypes]
        tmp1 = AnnData(
            np.concatenate(X, axis=0),
            obs=pd.DataFrame(celltypes_col, columns=["Celltype"]),
        )

        sc.pp.subsample(tmp1, fraction=0.5, random_state=42)
        sc.pp.normalize_total(tmp1, target_sum=1e6)
        sc.pp.log1p(tmp1)
        sc.tl.pca(tmp1)
        fig, ax = plt.subplots(nrows=1, ncols=1)
        sc.pl.pca(tmp1, color="Celltype", show=False, ax=ax)
        plt.savefig(
            os.path.join(
                self.simulation_folder,
                "scatterplot_pca_simulated_celltypes.pdf",
            ),
            bbox_inches="tight",
        )

    def simulate_per_batch(self, save=True):
        adata_orig = self.sc_adata.copy()
        adatas = []
        for batch in self.batches:
            self.sc_adata = adata_orig[
                adata_orig.obs[self.config["simulation_params"]["batch_col"]] == batch
            ]
            adatas.append(self.simulate(save=False))
        adata = adatas[0].concatenate(adatas[1:], batch_categories=self.batches)

        if save:
            adata.write(os.path.join(self.simulation_folder, "simulated.h5ad"))
            if self.config["simulation_params"]["generate_component_figures"]:
                self._generate_component_figures(adata)
        else:
            return adata

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
    # print(batch_col)
    # print(columns)
    if batch_col in columns:
        sim.batches = np.array(sim.sc_adata.obs[sim.config["simulation_params"]["batch_col"]].unique())
        if len(sim.batches) > 1:
            sim.simulate_per_batch(save=True)
        else:
            sim.simulate(save=True)
    else:
        print("Number of batches in single-cell data is 1. If this is incorrect, please specify name of the batch column as in the single-cell data object (.obs)")
        sim.simulate(save=True)
    sim.config["simulation_params"]["simulation_folder"] = os.path.join(sim.config["experiment_folder"], "simulation")
    if config["simulation_params"]["type"]=="bulk":
        sim.config["simulation_params"]["concentration"] = list(sim.config["simulation_params"]["concentration"])
    sim.config["deconv_params"]["reference"] = os.path.join(sim.config["simulation_params"]["simulation_folder"], "simulated.h5ad")
    save_dict_to_file(sim.config)
    
