import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc

import scipy


def convert_h5ad_to_df(adata):
    if type(adata.X) == np.ndarray:
        X = adata.X
    elif type(adata.X) == np.matrix:
        X = np.array(adata.X)
    elif (
        type(adata.X) == scipy.sparse.csr.csr_matrix
        or type(adata.X) == scipy.sparse.coo.coo_matrix
    ):
        X = np.array(adata.X.todense())
    else:
        sys.exit(
            "adata.X is of type {} which is not implemented.".format(type(adata.X))
        )
    genes = adata.var_names.tolist()
    df = pd.DataFrame(X, columns=genes, index=adata.obs["cell.type"].tolist())
    if "subject" in adata.obs.columns:
        subjects = pd.DataFrame(
            adata.obs.subject.tolist(), index=list(range(df.shape[0]))
        )
        return df, df_subjects
    else:
        return df


def save_to_csv(df, savepath, savename, sep, remove_duplicates=True, genes_in_row=True):
    if remove_duplicates:
        df = df.loc[:, ~df.columns.duplicated(keep="first")]
    celltypes = df.index.tolist()
    df.index = list(range(df.shape[0]))
    df_celltypes = pd.DataFrame(celltypes, index=df.index)
    if genes_in_row:
        df = df.T
    df.to_csv(os.path.join(savepath, savename + "_counts.txt"), sep)
    df_celltypes.to_csv(os.path.join(savepath, savename + "_celltypes.txt"), sep)


def save_sc(ref, savedir, method):
    savepath = os.path.join(savedir, method)
    if "/" in ref:
        savename = ref.split(".h5ad")[0].split("/")[-1]
    else:
        savename = ref.split(".h5ad")[0]
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    adata = sc.read(ref)
    if "subject" in adata.obs.columns:
        df, df_subjects = convert_h5ad_to_df(adata)
        df_subjects.to_csv(os.path.join(savepath, savename + "_subjects.txt"), sep=",")
    else:
        df = convert_h5ad_to_df(adata)
    if method == "CS_datasets":
        sep = "\t"
    else:
        sep = ","
    if "method" == "CS_datasets":
        savename = savename + "_cs"
    save_to_csv(df, savepath, savename, sep)


def save_test(test_path, test_format, savedir, method, remove_duplicates=True):
    if test_format == "txt":
        df = pd.read_table(test_path, index_col=0)
    if remove_duplicates:
        df = df.loc[~df.index.duplicated(keep="first")]
    if method == "CS_datasets":
        sep = "\t"
        df.index.name = "GeneSymbol"
    else:
        sep = ","
    savepath = os.path.join(savedir, method)
    if "/" in test_path:
        if "h5ad" in test_path:
            savename = test_path.split(".h5ad")[0].split("/")[-1]
        else:
            savename = test_path.split(".txt")[0].split("/")[-1]
    else:
        if "h5ad" in test_path:
            savename = test_path.split(".h5ad")[0]
        else:
            savename = test_path.split(".txt")[0]
    if "method" == "CS_datasets":
        savename = savename + "_cs"
    df.to_csv(os.path.join(savepath, savename + "_test.txt"), sep)
