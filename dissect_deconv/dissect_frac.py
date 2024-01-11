import os
import sys
import numpy as np
import pandas as pd
import tensorflow as tf
from dissect_deconv.utils.network_fn import network1 as network
from dissect_deconv.utils.network_fn import loss
from tqdm import tqdm
import random
from dissect_deconv.utils.utils_fn import normalize_per_batch, reproducibility, ccc_fn
from sklearn.metrics import mean_squared_error


def run_dissect_frac(config):
    dataset_path = os.path.join(config["experiment_folder"], "datasets")
    if not os.path.exists(dataset_path):
        sys.exit(
            "Path {} does not exist. Please run prepare_data.py before.".format(
                dataset_path
            )
        )

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

    j = 0
    for seed in range(len(config["deconv_params"]["models"])):
        print("Starting training model {}".format(j))
        # reproducibility(seed)

        minval, maxval = config["deconv_params"]["alpha_range"]
        # Create dataset iterators
        np.random.seed(seed)
        tf.random.set_seed(seed)
        dataset = tf.data.Dataset.from_tensor_slices((X_sim_np, y_sim_np, X_real_np))
        dataset = (
            dataset.shuffle(1000)
            .repeat()
            .batch(batch_size=config["deconv_params"]["network_params"]["batch_size"])
        )
        dataset_iter = iter(dataset)

        if config["deconv_params"]["network_params"]["hidden_activation"] == "relu6":
            config["deconv_params"]["network_params"]["hidden_activation"] = tf.nn.relu6
        else:
            config["deconv_params"]["network_params"]["hidden_activation"] = tf.nn.relu
        model = network(
            config["deconv_params"]["network_params"], n_celltypes, n_features, training=True
        )

        # Start training
        pbar = tqdm(range(config["deconv_params"]["network_params"]["n_steps"]))
        step = 0
        optimizer = tf.keras.optimizers.Adam(
            learning_rate=config["deconv_params"]["network_params"]["lr"]
        )
        rs, rmses, cccs, avgrs, avgrmses, avgcccs = [], [], [], [], [], []
        for i in pbar:
            X_sim, y_sim, X_real = dataset_iter.get_next()

            # seed = np.random.uniform()
            alpha = tf.random.uniform(
                [1], minval=minval, maxval=maxval, dtype=tf.dtypes.float32, name="alpha"
            )
            if config["deconv_params"]["mix"] == "rrm":
                X_real_s = tf.random.shuffle(X_real)
                X_mix = alpha * X_real + (1 - alpha) * X_real_s
            else:
                X_mix = alpha * X_real + (1 - alpha) * X_sim

            X_real, X_sim, X_mix = (
                normalize_per_batch(X_real, n_features),
                normalize_per_batch(X_sim, n_features),
                normalize_per_batch(X_mix, n_features),
            )
            if config["deconv_params"]["mix"] == "rrm":
                X_real_s = normalize_per_batch(X_real_s, n_features)
            reg_losses, cons_losses, total_losses = [], [], []

            if i == 0:
                y_hat_sim, y_hat_real, y_hat_mix = (
                    model(X_sim),
                    model(X_real),
                    model(X_mix),
                )
                print("Network architecture -")
                print(model.summary())

            with tf.GradientTape() as tape:
                if config["deconv_params"]["mix"] == "rrm":
                    y_hat_sim, y_hat_real, y_hat_mix, y_hat_real_s = (
                        model(X_sim, training=True),
                        model(X_real, training=True),
                        model(X_mix, training=True),
                        model(X_real_s, training=True),
                    )
                    reg_loss, cons_loss = loss(
                        config["deconv_params"]["network_params"]["loss"],
                        y_hat_sim,
                        y_sim,
                        y_hat_real,
                        y_hat_mix,
                        alpha,
                        y_hat_real_s,
                    )
                else:
                    y_hat_sim, y_hat_real, y_hat_mix = (
                        model(X_sim, training=True),
                        model(X_real, training=True),
                        model(X_mix, training=True),
                    )
                    reg_loss, cons_loss = loss(
                        config["deconv_params"]["network_params"]["loss"],
                        y_hat_sim,
                        y_sim,
                        y_hat_real,
                        y_hat_mix,
                        alpha,
                    )
                if step < 2000:
                    loss_ = reg_loss
                elif step in range(2000, 4000):
                    lambda_ = 15
                    if (
                        config["deconv_params"]["sig_matrix"]
                        or config["deconv_params"]["test_dataset_type"] == "spatial_sparse"
                    ):
                        lambda_ = 0.15
                    loss_ = reg_loss + lambda_ * cons_loss
                elif step >= 4000:
                    lambda_ = 10
                    if (
                        config["deconv_params"]["sig_matrix"]
                        or config["deconv_params"]["test_dataset_type"] == "spatial_sparse"
                    ):
                        lambda_ = 0.10
                    loss_ = reg_loss + lambda_ * cons_loss
            grads = tape.gradient(loss_, model.trainable_weights)
            optimizer.apply_gradients(zip(grads, model.trainable_weights))

            del tape

            total_losses.append(loss_)
            reg_losses.append(reg_loss)
            cons_losses.append(cons_loss)
            evaluate = False
            if evaluate:
                yt, yt_ = np.array(gt), model.predict(
                    normalize_per_batch(X_real_test, n_features)
                )
                yt, yt_ = gt, model.predict(
                    normalize_per_batch(X_real_test, n_features)
                )
                yt_ = pd.DataFrame(yt_, columns=celltypes)
                celltypes_select = [col for col in gt.columns if col != "Unknown"]
                yt, yt_ = np.array(gt[celltypes_select]), np.array(
                    yt_[celltypes_select]
                )
                s = (yt.shape[0] * yt.shape[1],)
                y, y_ = yt.reshape(s), yt_.reshape(s)

                r = np.corrcoef(y, y_)[0, 1]
                rmse = mean_squared_error(y, y_, squared=False)
                ccc = ccc_fn(y, y_)

                rs.append(r)
                rmses.append(rmse)
                cccs.append(ccc)

                avgr, avgrmse, avgccc = 0, 0, 0
                for k in range(yt.shape[1]):
                    if celltypes[k] != "Unknown":
                        avgr += np.corrcoef(yt[:, k], yt_[:, k])[0, 1]
                        avgrmse += mean_squared_error(
                            yt[:, k], yt_[:, k], squared=False
                        )
                        avgccc += ccc_fn(yt[:, k], yt_[:, k])
                avgr = avgr / (yt.shape[1] - 1)
                avgrmse = avgrmse / (yt.shape[1] - 1)
                avgccc = avgccc / (yt.shape[1] - 1)

                avgrs.append(avgr)
                avgrmses.append(avgrmse)
                avgcccs.append(avgccc)

                pbar.set_description(
                    "step: %d| loss: %.4f | r: %.4f | Avg r: %.4f"
                    % (step, loss_, r, avgr)
                )
            step += 1
            pbar.set_description(
                "step: %d| loss: %.4f"
                % (step, loss_)
            )

        model_path = os.path.join(config["experiment_folder"], "model_{}".format(j))
        model.save(model_path)

        model_p = tf.keras.models.load_model(model_path)

        print("Running deconvolution")
        y_hat = model_p.predict(normalize_per_batch(X_real_test, n_features))
        df_y_hat = pd.DataFrame(y_hat, columns=celltypes)
        df_y_hat.index = sample_names
        results_path = os.path.join(
            config["experiment_folder"], "dissect_fractions_{}.txt".format(j)
        )
        df_y_hat.to_csv(results_path, sep="\t")

        # metrics = ["r", "rmse", "ccc", "avgr", "avgrmse", "avgccc"]
        # metrics_vals = [rs, rmses, cccs, avgrs, avgrmses, avgcccs]
        # df_metrics = pd.DataFrame(columns=metrics,
        #                         index=range(config["deconv_params"]["network_params"]["n_steps"]))
        # for m in range(len(metrics)):
        #    df_metrics[metrics[m]] = metrics_vals[m]
        # df_metrics.to_csv(os.path.join(config["experiment_folder"], "per_step_metrics_{}.txt".format(j)),
        #                sep="\t")

        print("Estimated proportions are saved at {}.".format(results_path))
        if i < 5:
            print("Starting training model {}".format(i))
        j += 1


if __name__ == "__main__":
    from configs.main_config import config
    run_dissect_frac(config)
