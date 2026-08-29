import os
import sys
import numpy as np
import pandas as pd
import tensorflow as tf
from dissect.utils.network_fn import network1 as network
from dissect.utils.network_fn import loss
from tqdm import tqdm
import random
from dissect.utils.utils_fn import normalize_per_batch, reproducibility, ccc_fn
from sklearn.metrics import mean_squared_error
import logging


def _quiet_tf():
    """Silence TensorFlow's C++ and Python logging."""
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    logging.getLogger('tensorflow').setLevel(logging.FATAL)


def _dataset_path(config):
    """Return the prepared datasets folder, exiting if it does not exist."""
    dataset_path = os.path.join(config["experiment_folder"], "datasets")
    if not os.path.exists(dataset_path):
        sys.exit(
            "Path {} does not exist. Please run prepare_data.py before.".format(
                dataset_path
            )
        )
    return dataset_path


def load_datasets(config):
    """
    Read the prepared arrays and labels from the experiment's datasets folder.

    Returns a dict with the simulated training data (X_sim, y_sim), the real
    training and test expression matrices, the sample names, the cell types,
    and the derived feature and cell-type counts.
    """
    dataset_path = _dataset_path(config)

    print("Loading prepared datasets...")
    X_real_np = np.load(os.path.join(dataset_path, "X_real_train.npy"))
    X_sim_np = np.load(os.path.join(dataset_path, "X_sim.npy"))
    y_sim_np = np.load(os.path.join(dataset_path, "y_sim.npy"))
    X_real_test = np.load(os.path.join(dataset_path, "X_real_test.npy"))
    sample_names = pd.read_table(
        os.path.join(dataset_path, "sample_names.txt"), index_col=0
    ).index.tolist()
    celltypes = pd.read_table(
        os.path.join(dataset_path, "celltypes.txt"), index_col=0
    ).index.tolist()

    X_real_np, X_sim_np, y_sim_np = (
        np.array(X_real_np, dtype=np.float32),
        np.array(X_sim_np, dtype=np.float32),
        np.array(y_sim_np, dtype=np.float32),
    )
    X_real_test = np.array(X_real_test, dtype=np.float32)
    n_features = X_sim_np.shape[1]
    n_celltypes = len(celltypes)

    return {
        "X_real_np": X_real_np,
        "X_sim_np": X_sim_np,
        "y_sim_np": y_sim_np,
        "X_real_test": X_real_test,
        "sample_names": sample_names,
        "celltypes": celltypes,
        "n_features": n_features,
        "n_celltypes": n_celltypes,
    }


def _resolve_hidden_activation(config):
    """
    Replace the hidden_activation string in the config with the TensorFlow
    function it names: "relu6" becomes tf.nn.relu6, anything else tf.nn.relu.

    Does nothing once the entry is already a function, so the resolved value
    survives repeated calls on the same config object.
    """
    params = config["deconv_params"]["network_params"]
    if isinstance(params["hidden_activation"], str):
        if params["hidden_activation"] == "relu6":
            params["hidden_activation"] = tf.nn.relu6
        else:
            params["hidden_activation"] = tf.nn.relu


def model_paths(config):
    """Return the export paths of the ensemble models, in ensemble order."""
    return [
        os.path.join(config["experiment_folder"], "model_{}".format(j))
        for j in range(len(config["deconv_params"]["models"]))
    ]


def dissect_train(config, data=None):
    """
    Train the DISSECT ensemble and export each model to model_{j}
    inside the experiment folder.

    If use_consistency=True:
        Normal DISSECT training is performed.

        Steps 0-1999:
            regression loss only

        Steps 2000-3999:
            regression loss + 15 * consistency loss

        Steps >= 4000:
            regression loss + 10 * consistency loss

        For signature matrices / spatial_sparse:
            weights 0.15 and 0.10 are used instead.

    If use_consistency=False:
        A supervised baseline is trained exclusively on simulated
        samples and their known proportions.

        No real samples are used during training.
        No alpha values are generated.
        No online mixtures are generated.
        No real or mixed samples are normalized.
        No real/mixed forward passes are performed.
        No consistency loss is calculated.
    """

    _quiet_tf()

    if data is None:
        data = load_datasets(config)

    X_real_np = data["X_real_np"]
    X_sim_np = data["X_sim_np"]
    y_sim_np = data["y_sim_np"]
    X_real_test = data["X_real_test"]
    celltypes = data["celltypes"]
    n_features = data["n_features"]
    n_celltypes = data["n_celltypes"]

    network_params = config["deconv_params"]["network_params"]

    use_consistency = network_params.get(
        "use_consistency",
        True
    )

    # Resolve activation once before the ensemble loop
    _resolve_hidden_activation(config)

    print(
        "Consistency regularization: {}".format(
            "ON"
            if use_consistency
            else "OFF (supervised baseline)"
        )
    )

    if use_consistency:
        print(
            "Training mode: DISSECT simulation-real "
            "consistency training"
        )
    else:
        print(
            "Training mode: supervised simulated-data baseline"
        )

    saved_models = []

    j = 0

    # =========================================================
    # Train ensemble members
    # =========================================================
    for seed in range(
        len(config["deconv_params"]["models"])
    ):

        print(
            "Starting training model {}".format(j)
        )

        # -----------------------------------------------------
        # Reproducibility
        # -----------------------------------------------------
        np.random.seed(seed)
        tf.random.set_seed(seed)

        # -----------------------------------------------------
        # Dataset iterator
        #
        # Baseline does not even include X_real in the
        # TensorFlow training dataset.
        # -----------------------------------------------------
        if use_consistency:

            dataset = tf.data.Dataset.from_tensor_slices(
                (
                    X_sim_np,
                    y_sim_np,
                    X_real_np
                )
            )

        else:

            dataset = tf.data.Dataset.from_tensor_slices(
                (
                    X_sim_np,
                    y_sim_np
                )
            )

        dataset = (
            dataset
            .shuffle(1000)
            .repeat()
            .batch(
                batch_size=network_params[
                    "batch_size"
                ]
            )
        )

        dataset_iter = iter(dataset)

        # -----------------------------------------------------
        # Model
        # -----------------------------------------------------
        model = network(
            network_params,
            n_celltypes,
            n_features,
            training=True
        )

        # -----------------------------------------------------
        # Optimizer
        # -----------------------------------------------------
        optimizer = tf.keras.optimizers.Adam(
            learning_rate=network_params["lr"]
        )

        # -----------------------------------------------------
        # Training setup
        # -----------------------------------------------------
        pbar = tqdm(
            range(
                network_params["n_steps"]
            )
        )

        step = 0

        rs = []
        rmses = []
        cccs = []
        avgrs = []
        avgrmses = []
        avgcccs = []

        # Alpha is only relevant when consistency is enabled
        if use_consistency:

            minval, maxval = (
                config["deconv_params"]["alpha_range"]
            )

        # =====================================================
        # Training loop
        # =====================================================
        for i in pbar:

            # =================================================
            # NORMAL DISSECT
            # =================================================
            if use_consistency:

                X_sim, y_sim, X_real = (
                    dataset_iter.get_next()
                )

                # ---------------------------------------------
                # Generate alpha
                # ---------------------------------------------
                alpha = tf.random.uniform(
                    [1],
                    minval=minval,
                    maxval=maxval,
                    dtype=tf.dtypes.float32,
                    name="alpha"
                )

                # ---------------------------------------------
                # Generate online mixture
                # ---------------------------------------------
                if (
                    config["deconv_params"]["mix"]
                    == "rrm"
                ):

                    X_real_s = tf.random.shuffle(
                        X_real
                    )

                    X_mix = (
                        alpha * X_real
                        + (1 - alpha) * X_real_s
                    )

                else:

                    X_mix = (
                        alpha * X_real
                        + (1 - alpha) * X_sim
                    )

                # ---------------------------------------------
                # Normalize simulated, real and mixed data
                # ---------------------------------------------
                X_real = normalize_per_batch(
                    X_real,
                    n_features
                )

                X_sim = normalize_per_batch(
                    X_sim,
                    n_features
                )

                X_mix = normalize_per_batch(
                    X_mix,
                    n_features
                )

                if (
                    config["deconv_params"]["mix"]
                    == "rrm"
                ):

                    X_real_s = normalize_per_batch(
                        X_real_s,
                        n_features
                    )

                # ---------------------------------------------
                # Forward passes + loss
                # ---------------------------------------------
                with tf.GradientTape() as tape:

                    y_hat_sim = model(
                        X_sim,
                        training=True
                    )

                    y_hat_real = model(
                        X_real,
                        training=True
                    )

                    y_hat_mix = model(
                        X_mix,
                        training=True
                    )

                    if (
                        config["deconv_params"]["mix"]
                        == "rrm"
                    ):

                        y_hat_real_s = model(
                            X_real_s,
                            training=True
                        )

                        reg_loss, cons_loss = loss(
                            network_params["loss"],
                            y_hat_sim,
                            y_sim,
                            y_hat_real,
                            y_hat_mix,
                            alpha,
                            y_hat_real_s,
                            use_consistency=True
                        )

                    else:

                        reg_loss, cons_loss = loss(
                            network_params["loss"],
                            y_hat_sim,
                            y_sim,
                            y_hat_real,
                            y_hat_mix,
                            alpha,
                            use_consistency=True
                        )

                    # -----------------------------------------
                    # Original DISSECT consistency schedule
                    # -----------------------------------------
                    if step < 2000:

                        loss_ = reg_loss

                    elif step < 4000:

                        lambda_ = 15

                        if (
                            config["deconv_params"][
                                "sig_matrix"
                            ]
                            or
                            config["deconv_params"][
                                "test_dataset_type"
                            ]
                            == "spatial_sparse"
                        ):

                            lambda_ = 0.15

                        loss_ = (
                            reg_loss
                            + lambda_ * cons_loss
                        )

                    else:

                        lambda_ = 10

                        if (
                            config["deconv_params"][
                                "sig_matrix"
                            ]
                            or
                            config["deconv_params"][
                                "test_dataset_type"
                            ]
                            == "spatial_sparse"
                        ):

                            lambda_ = 0.10

                        loss_ = (
                            reg_loss
                            + lambda_ * cons_loss
                        )

            # =================================================
            # SUPERVISED BASELINE
            # =================================================
            else:

                # ---------------------------------------------
                # Only simulated samples and labels
                # ---------------------------------------------
                X_sim, y_sim = (
                    dataset_iter.get_next()
                )

                # ---------------------------------------------
                # Only normalize simulated data
                # ---------------------------------------------
                X_sim = normalize_per_batch(
                    X_sim,
                    n_features
                )

                # ---------------------------------------------
                # One forward pass only
                # ---------------------------------------------
                with tf.GradientTape() as tape:

                    y_hat_sim = model(
                        X_sim,
                        training=True
                    )

                    # Same loss() function as normal DISSECT,
                    # but consistency is disabled.
                    reg_loss, cons_loss = loss(
                        network_params["loss"],
                        y_hat_sim,
                        y_sim,
                        use_consistency=False
                    )

                    loss_ = reg_loss

            # =================================================
            # Backpropagation
            # =================================================
            grads = tape.gradient(
                loss_,
                model.trainable_weights
            )

            optimizer.apply_gradients(
                zip(
                    grads,
                    model.trainable_weights
                )
            )

            del tape

            step += 1

            pbar.set_description(
                "step: %d| loss: %.4f"
                % (
                    step,
                    float(loss_.numpy())
                )
            )

        # =====================================================
        # Save ensemble member
        # =====================================================
        model_path = os.path.join(
            config["experiment_folder"],
            "model_{}".format(j)
        )

        model.save(model_path)

        saved_models.append(model_path)

        j += 1

    # =========================================================
    # Save list of trained models
    # =========================================================
    with open(
        os.path.join(
            config["experiment_folder"],
            "models.txt"
        ),
        "w"
    ) as fh:

        for p in saved_models:
            fh.write(p + "\n")

    return saved_models


def dissect_predict(config, data=None):
    """
    Deconvolve the real test data with the exported ensemble.

    Each model_{j} is loaded from the experiment folder and applied to the
    normalised test matrix, giving per-model cell-type fractions and the
    corresponding pre-softmax scores. Both are written as
    dissect_fractions_{j}.txt and dissect_scores_{j}.txt, then averaged across
    the ensemble into dissect_fractions.txt and dissect_scores.txt.

    Loads the datasets itself when data is not supplied, so it can be run on
    its own after dissect_train.
    """
    _quiet_tf()

    if data is None:
        data = load_datasets(config)

    X_real_test = data["X_real_test"]
    sample_names = data["sample_names"]
    celltypes = data["celltypes"]
    n_features = data["n_features"]

    for j, model_path in enumerate(model_paths(config)):
        model_p = tf.keras.models.load_model(model_path)

        # print("Running deconvolution")
        y_hat = model_p.predict(normalize_per_batch(X_real_test, n_features))
        df_y_hat = pd.DataFrame(y_hat, columns=celltypes)
        df_y_hat.index = sample_names
        results_path = os.path.join(
            config["experiment_folder"], "dissect_fractions_{}.txt".format(j)
        )
        df_y_hat.to_csv(results_path, sep="\t")

        ## Output score
        last_layer = model_p.layers[-1]
        last_layer_config = last_layer.get_config()
        weights, biases = last_layer.get_weights()
        last_layer_config['activation'] = None
        new_last_layer = tf.keras.layers.Dense.from_config(last_layer_config)
        new_last_layer.build(input_shape=last_layer.input_shape)
        new_last_layer.set_weights([weights, biases])
        pre_softmax_model = tf.keras.models.Model(inputs=model_p.input, outputs=new_last_layer(model_p.layers[-2].output))
        scores = pre_softmax_model.predict(normalize_per_batch(X_real_test, n_features))

        df_y_hat_scores = pd.DataFrame(scores, columns=celltypes)
        df_y_hat_scores.index = sample_names
        scores_path = os.path.join(
            config["experiment_folder"], "dissect_scores_{}.txt".format(j)
        )
        df_y_hat_scores.to_csv(scores_path, sep="\t")

    # print("ensembling")
    experiment_path = os.path.join(config["experiment_folder"])

    i = 0
    for i in range(len(config["deconv_params"]["models"])):
        df_curr = pd.read_table(
            os.path.join(
                config["experiment_folder"], "dissect_fractions_{}.txt".format(i)
            ),
            index_col=0,
        )
        if i == 0:
            df_ens = df_curr
        else:
            df_ens = df_ens + df_curr

    df_ens = df_ens / len(config["deconv_params"]["models"])
    savepath = os.path.join(config["experiment_folder"], "dissect_fractions.txt")
    print("Predictions are saved to {}".format(savepath))
    df_ens.to_csv(savepath, sep="\t")

    for i in range(len(config["deconv_params"]["models"])):
        df_curr = pd.read_table(
            os.path.join(
                config["experiment_folder"], "dissect_scores_{}.txt".format(i)
            ),
            index_col=0,
        )
        if i == 0:
            df_ens = df_curr
        else:
            df_ens = df_ens + df_curr

    df_ens = df_ens / len(config["deconv_params"]["models"])
    savepath = os.path.join(config["experiment_folder"], "dissect_scores.txt")
    print("Scores are saved to {}".format(savepath))
    df_ens.to_csv(savepath, sep="\t")


def run_dissect_frac(config):
    """
    Run the full DISSECT fraction workflow.

    Loads the prepared datasets once, trains and exports the ensemble, then
    deconvolves the test data and writes the per-model and ensembled fraction
    and score tables.
    """
    _quiet_tf()

    data = load_datasets(config)

    dissect_train(config, data=data)
    dissect_predict(config, data=data)


if __name__ == "__main__":
    from configs.main_config import config
    run_dissect_frac(config)