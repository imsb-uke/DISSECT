import os
import random
import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.metrics import mean_squared_error as mse


def rmse_fn(x, y):
    return mse(x, y, squared=False)


def pcor(x, y):
    return np.corrcoef(x, y)[0][1]


def ccc_fn(y_true, y_pred):
    cor = np.corrcoef(y_true, y_pred)[0][1]
    mean_true = np.mean(y_true)
    mean_pred = np.mean(y_pred)
    var_true = np.var(y_true)
    var_pred = np.var(y_pred)

    sd_true = np.std(y_true)
    sd_pred = np.std(y_pred)

    numerator = 2 * cor * sd_true * sd_pred

    denominator = var_true + var_pred + (mean_true - mean_pred) ** 2

    return numerator / denominator


def evaluate(df1, df2, savepath):
    common = list(set(df1.columns) & set(df2.columns))
    common = [celltype for celltype in common if celltype != "Unknown"]
    df1, df2 = df1[common], df2[common]
    df_results = pd.DataFrame(
        columns=["r", "ccc", "rmse"], index=common + ["avg", "overall"]
    )
    avg_r, avg_ccc, avg_rmse = 0, 0, 0
    for celltype in common:
        x, y = df1[celltype], df2[celltype]
        r, ccc, rmse = pcor(x, y), ccc_fn(x, y), rmse_fn(x, y)
        avg_r, avg_ccc, avg_rmse = avg_r + r, avg_ccc + ccc, avg_rmse + rmse
        df_results.loc[celltype, ["r", "ccc", "rmse"]] = r, ccc, rmse
    avg_r, avg_ccc, avg_rmse = (
        avg_r / len(common),
        avg_ccc / len(common),
        avg_rmse / len(common),
    )
    df_results.loc["avg", ["r", "ccc", "rmse"]] = avg_r, avg_ccc, avg_rmse
    n = df1.shape[0] * df1.shape[1]
    x, y = np.array(df1).reshape((n,)), np.array(df2).reshape((n,))
    overall_r, overall_ccc, overall_rmse = pcor(x, y), ccc_fn(x, y), rmse_fn(x, y)
    df_results.loc["overall", ["r", "ccc", "rmse"]] = (
        overall_r,
        overall_ccc,
        overall_rmse,
    )

    return df_results


def log1p(x):
    ones = tf.ones_like(x, name="ones")
    x1 = tf.math.add(x, ones)

    numerator = tf.math.log(x1)
    denominator = tf.math.log(tf.constant(2, dtype=numerator.dtype))

    return numerator / denominator


def normalize_per_batch(x, n_features, epsilon=1e-8):
    x1 = log1p(x)
    min_, max_ = tf.reduce_min(x1, axis=1), tf.reduce_max(x1, axis=1)
    min1, max1 = tf.tile(
        tf.expand_dims(min_, axis=1), tf.constant([1, n_features])
    ), tf.tile(tf.expand_dims(max_, axis=1), tf.constant([1, n_features]))
    x_normed = (x1 - min1) / (
        max1 - min1 + epsilon
    )  # epsilon to avoid dividing by zero

    return x_normed


def set_seeds(seed=42):
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    tf.random.set_seed(seed)
    np.random.seed(seed)


def set_global_determinism(seed=42):
    set_seeds(seed=seed)

    os.environ["TF_DETERMINISTIC_OPS"] = str(1)
    os.environ["TF_CUDNN_DETERMINISTIC"] = str(1)
    # Read more here - https://www.tensorflow.org/api_docs/python/tf/config/threading/set_inter_op_parallelism_threads
    tf.config.threading.set_inter_op_parallelism_threads(1)
    tf.config.threading.set_intra_op_parallelism_threads(1)


def reproducibility(seed):
    set_seeds(seed)
    set_global_determinism(seed)
