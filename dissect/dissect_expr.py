from __future__ import absolute_import, division, print_function, unicode_literals
import os
import sys
import random
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
import anndata as ad
import decimal

from anndata import AnnData
import shutil
from tqdm import tqdm
import json
import matplotlib.pyplot as plt

import warnings
from tensorflow.keras import layers

import tensorflow as tf

import numpy as np

from keras.layers import Input, Dense, Lambda
from keras.layers.merge import concatenate as concat
from keras.models import Model
from keras import backend as K
from keras.datasets import mnist
from tensorflow.keras.utils import to_categorical
import scipy

warnings.filterwarnings("ignore")


def run_dissect_expr(config):
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    test_format = config["deconv_params"]["test_dataset_format"]
    real_path = config["deconv_params"]["test_dataset"]

    data = sc.read(config["deconv_params"]["reference"])
    data.var_names_make_unique()

    data.obs = data.obs[[col for col in data.obs.columns if col not in ["ds", "batch"]]]
    fractions = np.array(data.obs, dtype=np.float32)

    if test_format == "h5ad":
        real_data = sc.read(real_path)
    else:
        real_data = AnnData(pd.read_table(real_path, index_col=0).T)

    ens_path = os.path.join(config["experiment_folder"], "dissect_fractions_ens.txt")
    zero_path = os.path.join(config["experiment_folder"], "dissect_fractions_0.txt")

    if not os.path.exists(ens_path):
        if not os.path.exists(zero_path):
            sys.exit(
                "Could not find prediction proportions for the real data. Did you run dissect.py. They are in experiment_folder specified in the config file and have file name dissect_fraczions_[model_number].txt). "
            )
        else:
            real_fractions = pd.read_table(zero_path, index_col=0)
    else:
        real_fractions = pd.read_table(ens_path, index_col=0)

    real_fractions.index = real_data.obs_names
    real_data.obs = real_fractions

    real_data.var_names_make_unique()

    sc.pp.normalize_total(real_data, target_sum=1e6)
    sc.pp.log1p(real_data)

    sc.pp.normalize_total(data, target_sum=1e6)
    sc.pp.log1p(data)

    for layer in data.layers:
        tmp = AnnData(data.layers[layer])
        sc.pp.normalize_total(tmp, target_sum=1e6)
        sc.pp.log1p(tmp)
        data.layers[layer] = tmp.X

    X_real = real_data.to_df()
    X_real = X_real.loc[:, ~X_real.columns.duplicated(keep="first")]

    common_genes = list(set(X_real.columns.tolist()) & set(data.var_names.tolist()))
    print("Estimating expression of {} common genes".format(len(common_genes)))

    real_data = real_data[:, common_genes]
    data = data[:, common_genes]

    X_real = np.array(real_data.X)

    X_celltypes = {}
    X_sim = np.array(data.X)
    for layer in data.layers:
        if scipy.sparse.issparse(data.layers[layer]):
            data.layers[layer] = data.layers[layer].toarray()
        X_celltypes[layer] = np.array(data.layers[layer])
        df = pd.DataFrame(X_celltypes.keys(), columns=["Celltype"])

    df["Celltype"] = df.Celltype.astype("category")
    hot_encoding = np.array(pd.get_dummies(df.Celltype))

    mean_sum_sim = np.mean(np.sum(X_sim, axis=1))
    X_real_normalized = X_real / np.sum(X_real, axis=1, keepdims=True) * mean_sum_sim
    X_sim_normalized = X_sim.copy()

    scaling_method = config["deconv_params"]["network_params"]["expr_scaling"]
    if scaling_method == "p99":
        max_val = np.percentile(
            np.concatenate([X_sim_normalized.flatten(), X_real_normalized.flatten()]),
            99
        )
    else:
        max_val = max([X_sim_normalized.max(), X_real_normalized.max()])

    X_sim = X_sim_normalized / max_val
    X_real = X_real_normalized / max_val

    for celltype in X_celltypes.keys():
        X_celltypes[celltype] = np.concatenate(
            [
                X_celltypes[celltype],
                np.array(data.obs[celltype]).reshape((data.shape[0], 1)),
            ],
            1,
        )

    X_sim_input = np.repeat(X_sim, repeats=df.shape[0], axis=0)
    labels = np.concatenate([hot_encoding] * X_sim.shape[0], axis=0)
    X_sim_input = np.concatenate([X_sim_input, labels], axis=1)

    tmp = X_celltypes.copy()
    for celltype in tmp.keys():
        tmp[celltype] = pd.DataFrame(
            tmp[celltype],
            index=["s_{}".format(i) for i in range(tmp[celltype].shape[0])],
        )

    tmp_ls = []
    print("Preparing training datasets.")
    for i in range(X_celltypes[list(X_celltypes.keys())[0]].shape[0]):
        tmp_ls.append(
            np.concatenate(
                [
                    X_celltypes[celltype][i, :].reshape(
                        (1, X_celltypes[celltype].shape[1])
                    )
                    for celltype in X_celltypes.keys()
                ],
                axis=0,
            )
        )
    tmp = np.concatenate(tmp_ls, axis=0)

    hot_encoding = np.array(pd.get_dummies(df.Celltype))
    labels = np.concatenate([hot_encoding] * X_real.shape[0], axis=0)

    X_sim_gt = tmp.copy()
    X_sim_gt = X_sim_gt / max_val

    X_real_input = np.repeat(X_real, repeats=df.shape[0], axis=0)
    labels = np.concatenate([hot_encoding] * X_real.shape[0], axis=0)
    X_real_input = np.concatenate([X_real_input, labels], axis=1)

    labels = np.repeat(np.array(hot_encoding), repeats=X_real.shape[0], axis=0)

    class Sampling(layers.Layer):
        def call(self, inputs):
            z_mean, z_log_var = inputs
            batch = tf.shape(z_mean)[0]
            dim = tf.shape(z_mean)[1]
            epsilon = tf.keras.backend.random_normal(
                shape=(batch, dim), mean=0.0, stddev=1.0
            )
            return z_mean + tf.exp(0.5 * z_log_var) * epsilon

    class Encoder(layers.Layer):
        def __init__(
            self,
            latent_dim=32,
            intermediate_dim1=512,
            intermediate_dim2=512,
            intermediate_dim3=256,
            name="encoder",
            **kwargs
        ):
            super(Encoder, self).__init__(name=name, **kwargs)
            self.dense_proj2 = layers.Dense(intermediate_dim2, activation="relu")
            self.dense_proj3 = layers.Dense(intermediate_dim3, activation="relu")
            self.dense_mean = layers.Dense(latent_dim)
            self.dense_log_var = layers.Dense(latent_dim)
            self.sampling = Sampling()

        def call(self, inputs):
            x = self.dense_proj2(inputs)
            x = self.dense_proj3(x)
            z_mean = self.dense_mean(x)
            z_log_var = self.dense_log_var(x)
            z = self.sampling((z_mean, z_log_var))
            return z_mean, z_log_var, z

    class Decoder(layers.Layer):
        def __init__(
            self,
            original_dim,
            intermediate_dim1=128,
            intermediate_dim2=256,
            intermediate_dim3=512,
            name="decoder",
            **kwargs
        ):
            super(Decoder, self).__init__(name=name, **kwargs)
            self.dense_proj5 = layers.Dense(intermediate_dim2, activation="relu")
            self.dense_proj6 = layers.Dense(intermediate_dim3, activation="relu")
            self.dense_output = layers.Dense(
                original_dim - labels.shape[1], activation="relu"
            )

        def call(self, inputs, labels):
            x = tf.keras.layers.Concatenate(axis=-1)([inputs, labels])
            x = self.dense_proj5(x)
            x = self.dense_proj6(x)
            return self.dense_output(x)

    class VariationalAutoEncoder(tf.keras.Model):
        def __init__(
            self,
            original_dim,
            intermediate_dim1=200,
            intermediate_dim2=200,
            intermediate_dim3=200,
            latent_dim=128,
            name="autoencoder",
            **kwargs
        ):
            super(VariationalAutoEncoder, self).__init__(name=name, **kwargs)
            self.original_dim = original_dim
            self.encoder = Encoder(
                latent_dim=latent_dim,
                intermediate_dim1=intermediate_dim1,
                intermediate_dim2=intermediate_dim2,
                intermediate_dim3=intermediate_dim3,
            )
            self.decoder = Decoder(
                original_dim,
                intermediate_dim1=intermediate_dim3,
                intermediate_dim2=intermediate_dim2,
                intermediate_dim3=intermediate_dim1,
            )

        def call(self, inputs, labels):
            z_mean, z_log_var, z = self.encoder(inputs)
            reconstructed = self.decoder(z, labels)
            kl_loss = -0.5 * tf.reduce_mean(
                z_log_var - tf.square(z_mean) - tf.exp(z_log_var) + 1
            )
            return reconstructed, kl_loss

    fract = X_sim_gt[:, -1]

    fract_real = np.array(real_fractions, dtype=np.float32).reshape(-1, 1)
    fract_real = fract_real.reshape((fract_real.shape[0],))

    X_sim_gt = X_sim_gt[:, :-1]

    batch_size = 128

    data = tf.data.Dataset.from_tensor_slices(
        (X_sim_input, X_sim_gt, X_sim_input[:, X_sim_gt.shape[1]:], fract)
    )
    data = data.repeat().batch(batch_size=batch_size)
    data_iter = iter(data)

    data1 = tf.data.Dataset.from_tensor_slices(
        (X_real_input, X_real_input[:, X_sim_gt.shape[1]:], fract_real)
    )
    data1 = data1.repeat().batch(batch_size=batch_size)
    data_iter1 = iter(data1)

    labels = X_real_input[:, X_sim_gt.shape[1]:]
    if not config["deconv_params"]["network_params"]["n_steps_expr"]:
        total_steps = 5000 * df.shape[0]
    else:
        total_steps = config["deconv_params"]["network_params"]["n_steps_expr"]
    n_epochs = int((total_steps * 128) / X_sim_input[0].shape[0])

    original_dim = X_sim_input.shape[1]
    n_steps = int(X_sim_input.shape[0] / batch_size)

    vae = VariationalAutoEncoder(original_dim)
    mse_loss_fn = tf.keras.losses.MeanSquaredError(
        reduction=tf.keras.losses.Reduction.SUM
    )
    loss_metric = tf.keras.metrics.Sum()

    optimizer = tf.keras.optimizers.Adam(learning_rate=1e-3)

    pbar = tqdm(range(n_epochs))

    lambda_2 = 0.1
    beta_vae = 0.01

    i = 0
    for epoch in pbar:
        total_loss = 0
        total_loss_kl = 0
        total_loss_recon = 0

        for step in range(n_steps):
            x_train_sim, x_train_sim_gt, labels_sim, true_fractions_sim = data_iter.get_next()
            x_train_real, labels_real, true_fractions_real = data_iter1.get_next()

            beta = tf.random.uniform(shape=(), minval=0.1, maxval=0.9)

            batch_size_actual = tf.minimum(tf.shape(x_train_sim)[0], tf.shape(x_train_real)[0])
            x_train_sim_batch = x_train_sim[:batch_size_actual]
            x_train_real_batch = x_train_real[:batch_size_actual]
            labels_sim_batch = labels_sim[:batch_size_actual]
            labels_real_batch = labels_real[:batch_size_actual]
            true_fractions_sim_batch = true_fractions_sim[:batch_size_actual]
            true_fractions_real_batch = true_fractions_real[:batch_size_actual]
            x_train_sim_gt_batch = x_train_sim_gt[:batch_size_actual]

            x_train_mix = (
                beta * x_train_real_batch[:, :x_train_sim_gt_batch.shape[1]]
                + (1 - beta) * x_train_sim_batch[:, :x_train_sim_gt_batch.shape[1]]
            )
            labels_mix = labels_sim_batch
            true_fractions_mix = beta * true_fractions_real_batch + (1 - beta) * true_fractions_sim_batch

            with tf.GradientTape() as tape:
                reconstructed_sim, kld = vae(
                    x_train_sim_batch[:, :x_train_sim_gt_batch.shape[1]],
                    labels=labels_sim_batch
                )
                recon_loss = mse_loss_fn(x_train_sim_gt_batch, reconstructed_sim)

                reconstructed_mix, _ = vae(x_train_mix, labels=labels_mix)
                reconstructed_real, _ = vae(
                    x_train_real_batch[:, :x_train_sim_gt_batch.shape[1]],
                    labels=labels_real_batch
                )

                recon_mix_weighted = reconstructed_mix * tf.expand_dims(true_fractions_mix, 1)
                recon_real_weighted = reconstructed_real * tf.expand_dims(true_fractions_real_batch, 1)
                recon_sim_weighted = reconstructed_sim * tf.expand_dims(true_fractions_sim_batch, 1)
                expected_mix = beta * recon_real_weighted + (1 - beta) * recon_sim_weighted

                consistency_loss = mse_loss_fn(recon_mix_weighted, expected_mix)

                loss = recon_loss + beta_vae * kld + lambda_2 * consistency_loss

                total_loss_kl += kld
                total_loss_recon += recon_loss
                total_loss += loss

                grads = tape.gradient(loss, vae.trainable_weights)
                optimizer.apply_gradients(zip(grads, vae.trainable_weights))

            del tape

        i += 1
        pbar.set_description("Epoch: %d| loss %f " % (epoch, total_loss.numpy()))

    est, _ = vae(
        X_real_input[:, :X_sim_gt.shape[1]],
        labels=X_real_input[:, X_sim_gt.shape[1]:],
    )
    est = est * max_val

    est_adata = AnnData(est.numpy())
    est_adata.var_names = common_genes
    est_adata.X[est_adata.X < 0] = 0

    est_adata.obs["Celltype"] = np.array(df.Celltype.tolist() * X_real.shape[0])

    mapping = dict(zip(df.Celltype.astype(str).tolist(), df.Celltype.tolist()))

    est_adata.obs["Celltype"] = est_adata.obs["Celltype"].astype(str).astype("category")
    est_adata.obs["Celltype"].replace(mapping, inplace=True)

    if test_format == "h5ad":
        real_data = sc.read(real_path)
    elif test_format == "txt":
        X_real = pd.read_table(real_path, index_col=0).T
        real_data = AnnData(X_real)

    real_data.var_names_make_unique()
    real_data = real_data[:, common_genes]
    for celltype in est_adata.obs.Celltype.unique():
        real_data.layers[celltype] = est_adata[est_adata.obs["Celltype"] == celltype].X

    savename = os.path.join(config["experiment_folder"], "est_expression_layered.h5ad")
    real_data.write(savename)
    print(
        f"Estimated gene expression per cell type is saved at {savename} and is included in the layers attribute of the anndata object."
    )

    real_data.obs["sample"] = real_data.obs.index.tolist()
    for layer in real_data.layers:
        if layer == list(real_data.layers.keys())[0]:
            adata_est = real_data.copy()
            adata_est.X = adata_est.layers[layer]
            adata_est.obs["cell_type"] = layer
            del adata_est.layers
        else:
            tmp = real_data.copy()
            tmp.X = tmp.layers[layer]
            tmp.obs["cell_type"] = layer
            del tmp.layers
            adata_est = adata_est.concatenate(tmp)
    adata_est.obs = adata_est.obs[[col for col in adata_est.obs.columns if col != "batch"]]

    adata_est.layers["scaled_counts"] = sc.pp.scale(adata_est, copy=True, max_value=10).X
    savename = os.path.join(config["experiment_folder"], "est_expression.h5ad")
    adata_est.write(savename)
    print(
        f"Estimated gene expression per cell type for all samples is saved at {savename}, corresponding samples are in sample attribute of the anndata.obs."
    )


if __name__ == "__main__":
    from configs.main_config import config
    run_dissect_expr(config)