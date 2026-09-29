from sklearn.metrics import r2_score, mean_squared_error

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import adelie as ad

import numpy as np
import pyarrow.feather as feather
import pandas as pd
import argparse

np.random.seed(42)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Train a Gaussian Elastic Net model for quantitative resistance"
    )
    parser.add_argument("--data", type=str, required=True)
    parser.add_argument("--metadata", type=str, required=True)
    parser.add_argument("--output_prefix", type=str, required=True)
    parser.add_argument("--train_prop", type=float, default=0.5)
    parser.add_argument("--n_threads", type=int, default=1)
    parser.add_argument("--grouped", action="store_true", default=False)
    parser.add_argument("--max_iters", type=float, default=1e5)
    parser.add_argument("--tol", type=float, default=1e-7)
    parser.add_argument("--alpha", type=float, default=1)
    parser.add_argument(
        "--min_samples",
        type=int,
        default=8,
        help="Minimum number of samples per target to keep",
    )
    return parser.parse_args()


def read_feather_data(file_path):
    return feather.read_feather(file_path)


def read_metadata(file_path):
    metadata = pd.read_table(file_path)
    if "sample_name" not in metadata.columns:
        raise ValueError("Metadata file must contain a sample_name column")
    return metadata


def get_metadata_columns(metadata):
    return metadata.columns[metadata.columns != "sample_name"]


def merge_and_split_data(data, metadata, metadata_col, train_prop=0.5):
    metadata = metadata[["sample_name", metadata_col]]
    merged = pd.merge(data, metadata, on="sample_name", how="left")
    merged = merged.dropna(subset=[metadata_col])

    # convert to numeric target
    try:
        merged[metadata_col] = pd.to_numeric(merged[metadata_col])
    except:
        return None, None, None, None, None

    X = merged.drop(["sample_name", metadata_col], axis=1)
    y = merged[metadata_col].to_numpy()

    # random split
    idx = np.random.permutation(len(y))
    n_train = int(len(y) * train_prop)

    train_idx = idx[:n_train]
    test_idx = idx[n_train:]

    X_train = X.iloc[train_idx]
    y_train = y[train_idx]

    if train_prop == 1:
        return np.asarray(X_train), None, y_train, None, X.columns

    X_test = X.iloc[test_idx]
    y_test = y[test_idx]

    return (
        np.asarray(X_train, dtype=np.float64),
        np.asarray(X_test, dtype=np.float64),
        y_train,
        y_test,
        X.columns,
    )


def get_group_ids(column_names):
    group_ids = []
    current_group = None

    for i, col in enumerate(column_names):
        parts = col.split("_")
        if len(parts) < 3:
            raise ValueError(f"Bad column format: {col}")

        group = parts[1]
        if group != current_group:
            group_ids.append(i)
            current_group = group

    return np.array(group_ids, dtype=np.int32)


def train_adelie_model(
    X_train, y_train, n_threads=1, group_ids=None, max_iters=1e5, tol=1e-7, alpha=1
):
    model = ad.GroupElasticNet(
        solver="cv_grpnet",
        family="gaussian"
    )

    max_iters = int(max_iters)

    if group_ids is not None:
        model.fit(
            X_train.astype(np.float64),
            y_train.astype(np.float64),
            n_threads=n_threads,
            groups=group_ids,
            max_iters=max_iters,
            tol=tol,
            alpha=alpha,
        )
    else:
        model.fit(
            X_train.astype(np.float64),
            y_train.astype(np.float64),
            n_threads=n_threads,
            max_iters=max_iters,
            tol=tol,
            alpha=alpha,
        )

    return model


def main():
    args = parse_args()

    output_pdf = args.output_prefix + "_confusion_matrices.pdf"
    output_coef = args.output_prefix + "_coefficients.tsv"

    data = read_feather_data(args.data)
    metadata = read_metadata(args.metadata)

    metadata_columns = get_metadata_columns(metadata)

    all_model_features = None

    with PdfPages(output_pdf) as pdf:
        for metadata_col in metadata_columns:
            print(f"Processing: {metadata_col}")

            X_train, X_test, y_train, y_test, feature_names = merge_and_split_data(
                data,
                metadata,
                metadata_col,
                train_prop=args.train_prop,
            )

            if X_train is None:
                print("Skipping (non-numeric target)\n")
                continue

            if args.grouped:
                group_ids = get_group_ids(feature_names)
            else:
                group_ids = None

            try:
                model = train_adelie_model(
                    X_train,
                    y_train,
                    n_threads=args.n_threads,
                    group_ids=group_ids,
                    tol=args.tol,
                    max_iters=args.max_iters,
                    alpha=args.alpha,
                )
            except Exception as e:
                print(f"Failed: {e}")
                continue

            # ===== TEST PERFORMANCE =====
            if X_test is not None:
                y_pred = model.predict(X_test)

                r2 = r2_score(y_test, y_pred)
                rmse = np.sqrt(mean_squared_error(y_test, y_pred))

                print(f"Test R²: {r2:.3f}")
                print(f"Test RMSE: {rmse:.3f}")
            else:
                r2 = None
                rmse = None

            # ===== TRAIN PERFORMANCE =====
            y_train_pred = model.predict(X_train)

            train_r2 = r2_score(y_train, y_train_pred)
            train_rmse = np.sqrt(mean_squared_error(y_train, y_train_pred))

            print(f"Train R²: {train_r2:.3f}")
            print(f"Train RMSE: {train_rmse:.3f}\n")

            # ===== COEFFICIENTS =====
            coef = model.coef_.toarray().flatten()

            df_features = pd.DataFrame({
                "feature": feature_names,
                "coefficient": coef
            })

            df_features = df_features[df_features["coefficient"] != 0]

            df_features["metadata_category"] = metadata_col
            df_features["r2"] = r2 if r2 is not None else "NA"
            df_features["rmse"] = rmse if rmse is not None else "NA"
            df_features["train_r2"] = train_r2
            df_features["train_rmse"] = train_rmse

            if all_model_features is None:
                all_model_features = df_features
            else:
                all_model_features = pd.concat([all_model_features, df_features])

            # ===== PLOT =====
            if X_test is not None:
                plt.figure()
                plt.scatter(y_test, y_pred, alpha=0.5)
                plt.xlabel("True")
                plt.ylabel("Predicted")
                plt.title(f"{metadata_col}\nR²: {r2:.3f}, RMSE: {rmse:.3f}")
                plt.tight_layout()
                pdf.savefig()
                plt.close()

    if all_model_features is not None:
        all_model_features.to_csv(output_coef, sep="\t", index=False)
    else:
        pd.DataFrame().to_csv(output_coef, sep="\t", index=False)


if __name__ == "__main__":
    main()
