import tensorflow as tf


def network1(config, n_celltypes, n_features, training=False):
    """
    Creates network from config file (config), number of celltypes (n_celltypes)
    """
    model = tf.keras.Sequential()
    for l in range(config["n_hidden_layers"]):
        # model.add(tf.keras.layers.Dense(config["hidden_units"][l], "relu"))
        model.add(
            tf.keras.layers.Dense(
                config["hidden_units"][l],
                activation=config["hidden_activation"],
                kernel_regularizer=tf.keras.regularizers.l1_l2(l1=1e-5, l2=1e-5),
            )
        )
        if config["dropout"]:
            model.add(tf.keras.layers.Dropout(config["dropout"][l]))

    model.add(
        tf.keras.layers.Dense(n_celltypes, activation=config["output_activation"])
    )

    return model


def network2(config, n_celltypes, n_features, training=False):
    """
    Creates network from config file (config), number of celltypes (n_celltypes)
    """
    # model = tf.keras.Sequential()
    input = tf.keras.layers.Input(shape=(n_features,))
    for l in range(config["n_hidden_layers"]):
        if l == 0:
            # x = tf.keras.layers.BatchNormalization()(input)
            x = tf.keras.layers.Dense(
                config["hidden_units"][l], activation=config["hidden_activation"]
            )(input)
        else:
            # y = tf.keras.layers.BatchNormalization()(x)
            # y = x
            x = tf.keras.layers.Dense(config["hidden_units"][l], activation=None)(x)
            y = tf.keras.layers.Activation(config["hidden_activation"])(x)
            y = tf.keras.layers.Dropout(0.2)(y, training=training)
            y = tf.keras.layers.Dense(config["hidden_units"][l], activation=None)(y)
            y = tf.keras.layers.Dropout(0.2)(y, training=training)
            x = tf.keras.layers.Add()([x, y])

    # x = tf.keras.layers.BatchNormalization()(x)
    x = tf.keras.layers.Activation(config["hidden_activation"])(x)
    x = tf.keras.layers.Dropout(0.2)(x, training=training)
    x = tf.keras.layers.Dense(n_celltypes, activation=config["output_activation"])(x)

    model = tf.keras.Model(inputs=input, outputs=x)
    return model


def loss(
    loss_fn,
    y_hat_sim,
    y_sim,
    y_hat_real=None,
    y_hat_mix=None,
    alpha=None,
    y_hat_real_s=None,
    use_consistency=True,
):
    """
    Returns regression and consistency losses.

    Parameters
    ----------
    loss_fn : str
        One of "kldivergence", "l2", or "l1".

    y_hat_sim : tf.Tensor
        Predictions for simulated samples.

    y_sim : tf.Tensor
        Ground-truth proportions for simulated samples.

    y_hat_real : tf.Tensor or None
        Predictions for real samples. Only required when
        use_consistency=True.

    y_hat_mix : tf.Tensor or None
        Predictions for mixed samples. Only required when
        use_consistency=True.

    alpha : tf.Tensor or None
        Mixing coefficient. Only required when
        use_consistency=True.

    y_hat_real_s : tf.Tensor or None
        Predictions for shuffled real samples when using RRM.

    use_consistency : bool
        If False, calculate only the supervised regression loss.
        If True, calculate both regression and consistency losses.
    """

    # ---------------------------------------------------------
    # Supervised regression loss
    # ---------------------------------------------------------
    if loss_fn == "kldivergence":

        KL = tf.keras.losses.KLDivergence()

        reg_loss = KL(
            y_sim,
            y_hat_sim
        )

    elif loss_fn == "l2":

        reg_loss = tf.reduce_mean(
            tf.math.square(
                y_sim - y_hat_sim
            )
        )

    elif loss_fn == "l1":

        reg_loss = tf.reduce_mean(
            tf.math.abs(
                y_sim - y_hat_sim
            )
        )

    else:

        raise ValueError(
            "loss_fn must be one of: "
            "'kldivergence', 'l2', 'l1'"
        )

    # ---------------------------------------------------------
    # Supervised-only baseline
    # ---------------------------------------------------------
    if not use_consistency:

        cons_loss = tf.constant(
            0.0,
            dtype=reg_loss.dtype
        )

        return reg_loss, cons_loss


    if y_hat_real is None:
        raise ValueError(
            "y_hat_real is required when use_consistency=True"
        )

    if y_hat_mix is None:
        raise ValueError(
            "y_hat_mix is required when use_consistency=True"
        )

    if alpha is None:
        raise ValueError(
            "alpha is required when use_consistency=True"
        )

    # ---------------------------------------------------------
    # Construct prediction-space mixup target
    # ---------------------------------------------------------
    if y_hat_real_s is not None:

        # Real-real mixture
        y_mix = (
            alpha * y_hat_real
            + (1 - alpha) * y_hat_real_s
        )

    else:

        # Simulation-real mixture
        y_mix = (
            alpha * y_hat_real
            + (1 - alpha) * y_hat_sim
        )

    # ---------------------------------------------------------
    # Consistency loss
    # ---------------------------------------------------------
    if loss_fn == "l1":

        cons_loss = tf.reduce_mean(
            tf.math.abs(
                y_mix - y_hat_mix
            )
        )

    else:


        cons_loss = tf.reduce_mean(
            tf.math.square(
                y_mix - y_hat_mix
            )
        )

    return reg_loss, cons_loss
