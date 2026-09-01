from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import OneHotEncoder

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import adelie as ad

import numpy as np
import scipy.stats as st

import pyarrow.feather as feather
import pandas as pd

import argparse
from pathlib import Path

np.random.seed(42)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Train a model to predict antibiotic resistance"
    )
    # Keep original separate train/test file structure
    parser.add_argument(
        "--train_features",
        type=str,
        help="Path to the training features file",
        required=True,
    )
    parser.add_argument(
        "--train_metadata",
        type=str,
        help="Path to the training metadata file",
        required=True,
    )
    parser.add_argument(
        "--test_features",
        type=str,
        help="Path to the test features file",
        required=True,
    )
    parser.add_argument(
        "--test_metadata",
        type=str,
        help="Path to the test metadata file",
        required=True,
    )
    parser.add_argument(
        "--output_prefix", type=str, help="Prefix for the output files", required=True
    )
    parser.add_argument(
        "--min_samples",
        type=int,
        default=28,
        help="Minimum number of samples per category to keep",
    )
    parser.add_argument(
        "--n_threads",
        type=int,
        default=1,
        help="Number of threads to use for training the model",
    )

    # New parameters from the generic version
    parser.add_argument(
        "--even_samples",
        action="store_true",
        help="Keep the same number of samples per class in the training set",
    )
    parser.add_argument(
        "--grouped",
        action="store_true",
        default=False,
        help="Use grouped elastic net based on feature name prefixes",
    )
    parser.add_argument(
        "--max_iters",
        type=float,
        default=1e5,
        help="Maximum number of iterations for the Adelie model training",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-7,
        help="Tolerance for the Adelie model training convergence",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=1,
        help="Alpha parameter for elastic net regularization (0 = ridge, 1 = lasso)",
    )
    return parser.parse_args()


def read_feather_data(file_path):
    return feather.read_feather(file_path)


def get_metadata_delimiter(file_path):
    suffix = Path(file_path).suffix.lower()
    if suffix == ".csv":
        return ","
    return "\t"


def read_metadata(file_path):
    # Read all metadata columns as strings to avoid dtype mismatches later
    metadata = pd.read_csv(file_path, sep=get_metadata_delimiter(file_path), dtype=str)
    if "sample_name" not in metadata.columns:
        raise ValueError("Metadata file must contain a sample_name column")
    # Ensure all columns are strings (defensive)
    metadata = metadata.astype(str)
    return metadata


def get_metadata_columns(metadata, min_samples=50):
    """
    Returns the columns of the metadata file except for the sample_name column
    Filters the columns so it will only return metadata that have more than two
    discrete values with greater than min_samples per category
    """
    filtered_metadata = metadata.loc[:, metadata.columns != "sample_name"]
    if min_samples > 0:
        filtered_metadata = filtered_metadata.loc[
            :, filtered_metadata.apply(lambda x: len(x.unique()) >= 2, axis=0)
        ]
        filtered_metadata = filtered_metadata.loc[
            :,
            filtered_metadata.apply(
                lambda x: sum(x.value_counts() >= min_samples) > 1, axis=0
            ),
        ]
    else:
        filtered_metadata = filtered_metadata.loc[
            :, filtered_metadata.apply(lambda x: len(x.unique()) >= 1, axis=0)
        ]
        filtered_metadata = filtered_metadata.loc[
            :,
            filtered_metadata.apply(
                lambda x: sum(x.value_counts() >= min_samples) >= 1, axis=0
            ),
        ]
    return filtered_metadata.columns


def merge_data(data, metadata, metadata_col, min_samples=50, even_samples=False):
    metadata = metadata[["sample_name", metadata_col]]
    metadata[metadata_col] = metadata[metadata_col].replace("nan", pd.NA)

    merged_data = pd.merge(data, metadata, on="sample_name", how="left")
    merged_data = merged_data.dropna(subset=[metadata_col])

    class_counts = merged_data[metadata_col].value_counts()
    class_counts = class_counts[class_counts >= min_samples]
    classes_to_keep = class_counts.index
    classes_to_keep = classes_to_keep[~pd.isna(classes_to_keep)]
    classes_to_keep = classes_to_keep[classes_to_keep != "nan"]

    if len(classes_to_keep) == 0:
        return None, None, None
    if len(classes_to_keep) < 2 and min_samples != 0:
        print("This logic is true")
        return None, None, None

    merged_data = merged_data[merged_data[metadata_col].isin(classes_to_keep)]

    if even_samples:
        num_to_keep = class_counts.min()
        indices_to_keep = (
            merged_data.groupby(metadata_col)
            .apply(
                lambda x: x.sample(n=num_to_keep, replace=False).index,
                include_groups=False,
            )
            .explode()
        )
        merged_data = merged_data.loc[indices_to_keep]

    X = merged_data.drop(["sample_name", metadata_col], axis=1)
    y = merged_data[metadata_col].to_numpy()
    return np.asfortranarray(X), y, X.columns


def get_group_ids(column_names):
    """
    Given a list of the column names for X, return a list of the starting
    index of each group based on the number following the first underscore.
    The column names are expected to be in the format [cluster|kmer]_<group>_<feature>_NUM
    """
    group_ids = []
    current_group = None
    for i, col in enumerate(column_names):
        parts = col.split("_")
        if len(parts) < 3:
            raise ValueError(
                f"Column name {col} does not have the expected format [cluster|kmer]_<group>_<feature>_..."
            )
        group = parts[1]
        if group != current_group:
            group_ids.append(i)
            current_group = group

    return np.array(group_ids, dtype=np.int32)


def train_adelie_model(
    X_train, y_train, n_threads=1, group_ids=None, max_iters=1e5, tol=1e-7, alpha=0.5
):
    oh = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
    y_train2 = oh.fit_transform(y_train[:, np.newaxis])

    X_train_wrap = ad.matrix.dense(
        np.asarray(X_train, dtype=np.float64), method="naive", n_threads=n_threads
    )

    max_iters = int(max_iters)

    model = ad.GroupElasticNet(solver="cv_grpnet", family="multinomial")
    if group_ids is not None:
        model.fit(
            X_train_wrap,
            y_train2.astype(np.float64),
            n_threads=n_threads,
            groups=group_ids,
            max_iters=max_iters,
            tol=tol,
            alpha=alpha,
        )
    else:
        model.fit(
            X_train_wrap,
            y_train2.astype(np.float64),
            n_threads=n_threads,
            max_iters=max_iters,
            tol=tol,
            alpha=alpha,
        )

    return model, oh


def append_confusion_matrix_rows(rows, metadata_col, matrix_name, cm, labels):
    labels = list(map(str, labels))
    for i, true_label in enumerate(labels):
        for j, predicted_label in enumerate(labels):
            count = int(cm[i, j]) if i < cm.shape[0] and j < cm.shape[1] else 0
            rows.append(
                {
                    "metadata_category": metadata_col,
                    "matrix": matrix_name,
                    "true_label": true_label,
                    "predicted_label": predicted_label,
                    "n_samples": count,
                }
            )


def main():
    args = parse_args()
    output_prefix = args.output_prefix
    output_pdf = output_prefix + "_confusion_matrices.pdf"
    output_coef = output_prefix + "_nonzero_coefficients.tsv"
    output_confusion_tsv = output_prefix + "_confusion_matrices.tsv"

    train_features = read_feather_data(args.train_features)
    train_metadata = read_metadata(args.train_metadata)
    test_features = read_feather_data(args.test_features)
    test_metadata = read_metadata(args.test_metadata)

    train_metadata_columns = get_metadata_columns(
        train_metadata, min_samples=args.min_samples
    )
    test_metadata_columns = get_metadata_columns(test_metadata, min_samples=0)
    metadata_columns = [
        col for col in train_metadata_columns if col in test_metadata_columns
    ]

    all_model_features = None
    confusion_matrix_rows = []

    with PdfPages(output_pdf) as pdf:
        for metadata_col in metadata_columns:
            print(f"Processing metadata column: {metadata_col}")
            print()

            X_train, y_train, model_features = merge_data(
                train_features,
                train_metadata,
                metadata_col,
                min_samples=args.min_samples,
                even_samples=args.even_samples,
            )
            X_test, y_test, _ = merge_data(
                test_features, test_metadata, metadata_col, min_samples=0
            )

            if X_test is not None:
                test_classes_to_keep = np.isin(y_test, np.unique(y_train))
                X_test = X_test[test_classes_to_keep]
                y_test = y_test[test_classes_to_keep]

            if X_train is None or X_test is None:
                print(
                    f"Skipping {metadata_col} as there are not enough samples after merging and filtering..."
                )
                print()
                continue

            num_classes = len(np.unique(y_train))
            print(f"Number of classes for {metadata_col}: {num_classes}")

            # Set group ids based on feature names if --grouped is supplied
            if args.grouped and num_classes < 4:
                group_ids = get_group_ids(model_features)
                print(f"Using grouped elastic net with {len(group_ids)} groups.")
            else:
                print("Not using grouped elastic net.")
                if args.grouped:
                    print(
                        f"Skipping grouped elastic net for {metadata_col} as it has {num_classes} classes (must be less than 4)."
                    )
                group_ids = None

            try:
                model, oh = train_adelie_model(
                    X_train,
                    y_train,
                    n_threads=args.n_threads,
                    group_ids=group_ids,
                    tol=args.tol,
                    max_iters=args.max_iters,
                    alpha=args.alpha,
                )
            except Exception as e:
                print(f"Failed to train model for {metadata_col}: {e}")
                continue

            # Test predictions
            print(X_test)
            yhat = model.predict(X_test.astype(np.float64))
            if len(np.unique(yhat)) < 2:
                print(f"Test predictions for {metadata_col} are all of one class.")
                unique_class = np.unique(yhat)[0]
                yhat_2d = np.zeros((y_test.size, len(oh.categories_[0])))
                yhat_2d[:, unique_class] = 1
                y_pred = oh.inverse_transform(yhat_2d).flatten()
            else:
                yhat_2d = np.zeros((yhat.size, len(oh.categories_[0])))
                yhat_2d[np.arange(yhat.size), yhat] = 1
                yhat = yhat_2d
                y_pred = oh.inverse_transform(yhat).flatten()

            cm = confusion_matrix(y_test, y_pred, labels=oh.categories_[0])
            print(f"Test confusion matrix for {metadata_col}")
            print(cm)

            # Train predictions for accuracy tracking
            yhat_train = model.predict(X_train.astype(np.float64))
            if len(np.unique(yhat_train)) < 2:
                print(f"Train predictions for {metadata_col} are all of one class.")
                unique_class_train = np.unique(yhat_train)[0]
                yhat_train_2d = np.zeros((y_train.size, len(oh.categories_[0])))
                yhat_train_2d[:, unique_class_train] = 1
                y_train_pred = oh.inverse_transform(yhat_train_2d).flatten()
            else:
                yhat_train_2d = np.zeros((yhat_train.size, len(oh.categories_[0])))
                yhat_train_2d[np.arange(yhat_train.size), yhat_train] = 1
                yhat_train = yhat_train_2d
                y_train_pred = oh.inverse_transform(yhat_train).flatten()

            cm_train = confusion_matrix(y_train, y_train_pred, labels=oh.categories_[0])
            print(f"Train confusion matrix for {metadata_col}")
            print(cm_train)

            append_confusion_matrix_rows(
                confusion_matrix_rows,
                metadata_col,
                "test",
                cm,
                oh.categories_[0],
            )
            append_confusion_matrix_rows(
                confusion_matrix_rows,
                metadata_col,
                "train",
                cm_train,
                oh.categories_[0],
            )

            # Extract coefficients
            coef = model.coef_
            metadata_categories = oh.categories_[0]
            model_features_expanded = [
                f"{feature}+{category}"
                for feature in model_features
                for category in metadata_categories
            ]

            model_features_df = pd.DataFrame(
                model_features_expanded, columns=["feature"]
            )
            model_features_df["coefficient"] = coef.toarray().flatten()
            model_features_df = model_features_df[model_features_df["coefficient"] != 0]

            model_features_df["feature"] = model_features_df["feature"].str.split("+")
            model_features_df["category"] = model_features_df["feature"].str[1]
            model_features_df["feature"] = model_features_df["feature"].str[0]

            model_features_df = (
                model_features_df.groupby("feature")
                .agg({"coefficient": list, "category": list})
                .reset_index()
            )
            model_features_df["classes"] = model_features_df["category"].apply(
                lambda x: "[" + ",".join(map(str, x)) + "]"
            )
            model_features_df["coefficients"] = model_features_df["coefficient"].apply(
                lambda x: "[" + ",".join(map(str, x)) + "]"
            )

            model_features_df["metadata_category"] = metadata_col

            # Calculate test metrics
            if cm.shape[0] > 2:
                accuracy = np.trace(cm) / np.sum(cm)
                specificity = None
                sensitivity = None
            else:
                accuracy = (cm[0][0] + cm[1][1]) / sum(sum(row) for row in cm)
                specificity = cm[0][0] / (cm[0][0] + cm[0][1])
                sensitivity = cm[1][1] / (cm[1][0] + cm[1][1])

            # Calculate train metrics
            train_accuracy = np.trace(cm_train) / np.sum(cm_train)

            if specificity is None:
                print(f"Test Accuracy: {accuracy:.2f}")
            else:
                print(
                    f"Test Accuracy: {accuracy:.2f}, Specificity: {specificity:.2f}, Sensitivity: {sensitivity:.2f}"
                )

            print(f"Train Accuracy: {train_accuracy:.2f}")

            model_features_df["accuracy"] = accuracy
            model_features_df["train_accuracy"] = train_accuracy
            model_features_df["specificity"] = (
                specificity if specificity is not None else "NA"
            )
            model_features_df["sensitivity"] = (
                sensitivity if sensitivity is not None else "NA"
            )

            # Format confusion matrix for output
            out_cm = [map(str, row) for row in cm]
            out_cm = "{" + ";".join([",".join(row) for row in out_cm]) + "}"
            model_features_df["confusion_matrix"] = out_cm

            model_features_df = model_features_df[
                [
                    "metadata_category",
                    "feature",
                    "accuracy",
                    "train_accuracy",
                    "sensitivity",
                    "specificity",
                    "confusion_matrix",
                    "classes",
                    "coefficients",
                ]
            ]

            if all_model_features is None:
                all_model_features = model_features_df
            else:
                all_model_features = pd.concat(
                    [all_model_features, model_features_df], axis=0
                )

            # Plot confusion matrix
            plt.figure()
            plt.imshow(
                cm / cm.sum(axis=1)[:, np.newaxis], cmap="viridis", vmin=0, vmax=1
            )
            plt.colorbar()
            for i in range(cm.shape[0]):
                for j in range(cm.shape[1]):
                    plt.text(j, i, f"{cm[i, j]}", ha="center", va="center")
            if cm.shape[0] > 2:
                plt.title(
                    f"{metadata_col}\nTest Accuracy: {accuracy:.2f}\nTrain Accuracy: {train_accuracy:.2f}"
                )
            else:
                plt.title(
                    f"{metadata_col}\nTest Accuracy: {accuracy:.2f}, Specificity: {specificity:.2f}, Sensitivity: {sensitivity:.2f}\nTrain Accuracy: {train_accuracy:.2f}"
                )
            plt.xlabel("Predicted")
            plt.ylabel("True")
            plt.xticks(range(cm.shape[1]), metadata_categories, rotation=45)
            plt.yticks(range(cm.shape[0]), metadata_categories, rotation=45)
            plt.tight_layout()
            pdf.savefig()
            plt.close()

            print()

        if all_model_features is None:
            plt.figure()
            plt.text(
                0.5, 0.5, "No metadata columns were processed", ha="center", va="center"
            )
            plt.axis("off")
            pdf.savefig()
            plt.close()

            columns = [
                "metadata_category",
                "feature",
                "accuracy",
                "train_accuracy",
                "sensitivity",
                "specificity",
                "confusion_matrix",
                "classes",
                "coefficients",
            ]
            pd.DataFrame(columns=columns).to_csv(output_coef, sep="\t", index=False)
        else:
            all_model_features.to_csv(
                output_coef, sep="\t", index=False, float_format="%.4f"
            )

    pd.DataFrame(
        confusion_matrix_rows,
        columns=[
            "metadata_category",
            "matrix",
            "true_label",
            "predicted_label",
            "n_samples",
        ],
    ).to_csv(output_confusion_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
