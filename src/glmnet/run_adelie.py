# from sklearn.model_selection import train_test_split
# from sklearn.metrics import r2_score
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import OneHotEncoder

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import adelie as ad

import numpy as np

# import scipy.stats as st

import pyarrow.feather as feather
import pandas as pd
from math import floor

# from os.path import basename
import argparse
from pathlib import Path

np.random.seed(42)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Train a model to predict antibiotic resistance"
    )
    parser.add_argument("--data", type=str, help="Path to the data file", required=True)
    parser.add_argument(
        "--metadata", type=str, help="Path to the metadata file", required=True
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
    parser.add_argument(
        "--balanced_test",
        action="store_true",
        help="Keep the same number of samples per class in the test set",
    )
    parser.add_argument(
        "--train_prop",
        type=float,
        default=0.5,
        help="Proportion of the data to use for training."
        "Grabs this proportion from the smallest class and then evenly samples that number from all other classes.",
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
    metadata = pd.read_csv(file_path, sep=get_metadata_delimiter(file_path))
    if "sample_name" not in metadata.columns:
        raise ValueError("Metadata file must contain a sample_name column")
    # mutate all columns to strings for categorical analysis
    metadata = metadata.apply(lambda x: x.astype(str))
    return metadata


def get_metadata_columns(metadata, min_samples=50):
    """
    Returns the columns of the metadata file except for the sample_name column
    Filters the columns so it will only return metadata that have more than two
    discrete values with greater than min_samples per category
    """
    filtered_metadata = metadata.loc[:, metadata.columns != "sample_name"]
    # filter out columns with less than 2 unique values
    filtered_metadata = filtered_metadata.loc[
        :, filtered_metadata.apply(lambda x: len(x.unique()) >= 2, axis=0)
    ]
    # only grab columns with two or more categories that have more than min_samples
    filtered_metadata = filtered_metadata.loc[
        :,
        filtered_metadata.apply(
            lambda x: sum(x.value_counts() > min_samples) > 1, axis=0
        ),
    ]
    return filtered_metadata.columns


def merge_and_split_data(
    data, metadata, metadata_col, min_samples=50, train_prop=0.5, balanced_test=False
):
    metadata = metadata[["sample_name", metadata_col]]
    merged_data = pd.merge(data, metadata, on="sample_name", how="left")
    merged_data = merged_data.dropna(subset=[metadata_col])

    # Check the distribution of classes for this metadata category
    class_counts = merged_data[metadata_col].value_counts()

    # Drop any classes with less than min_samples
    class_counts = class_counts[class_counts >= min_samples]
    classes_to_keep = class_counts.index
    classes_to_keep = classes_to_keep[~pd.isna(classes_to_keep)]
    classes_to_keep = classes_to_keep[classes_to_keep != "nan"]
    # if there are not >= 2 classes with >= min_samples samples, return None
    if len(classes_to_keep) < 2:
        return None, None, None, None, None
    merged_data = merged_data[merged_data[metadata_col].isin(classes_to_keep)]

    # Get the minimum number of samples per class
    # keep exactly half of the samples for each class for the training set
    # and keep the rest of the samples for the test set
    num_to_keep = class_counts.min()
    if pd.isna(num_to_keep):
        return None, None, None, None, None
    num_to_keep = floor(num_to_keep * train_prop)
    indices_to_keep = (
        merged_data.groupby(metadata_col)
        .apply(
            lambda x: x.sample(n=num_to_keep, replace=False).index, include_groups=False
        )
        .explode()
    )

    if train_prop == 1:
        # Split the data into training and test sets
        X_train = merged_data.drop(["sample_name", metadata_col], axis=1).loc[
            indices_to_keep
        ]
        model_features = X_train.columns
        y_train = merged_data[metadata_col].loc[indices_to_keep].to_numpy()

        return np.asfortranarray(X_train), None, y_train, None, model_features

    # If we want a balanced test set, keep the same number of samples per class in the test set
    if balanced_test:
        X_train = merged_data.drop(["sample_name", metadata_col], axis=1).loc[
            indices_to_keep
        ]
        model_features = X_train.columns
        y_train = merged_data[metadata_col].loc[indices_to_keep].to_numpy()

        test_indices = (
            merged_data.drop(indices_to_keep)
            .groupby(metadata_col)
            .apply(
                lambda x: x.sample(n=num_to_keep, replace=False).index,
                include_groups=False,
            )
            .explode()
        )
        X_test = merged_data.drop(["sample_name", metadata_col], axis=1).loc[
            test_indices
        ]
        y_test = merged_data[metadata_col].loc[test_indices].to_numpy()
    else:
        # Split the data into training and test sets
        X_train = merged_data.drop(["sample_name", metadata_col], axis=1).loc[
            indices_to_keep
        ]
        model_features = X_train.columns
        y_train = merged_data[metadata_col].loc[indices_to_keep].to_numpy()

        X_test = merged_data.drop(["sample_name", metadata_col], axis=1).drop(
            indices_to_keep
        )
        y_test = merged_data[metadata_col].drop(indices_to_keep).to_numpy()

    return (
        np.asfortranarray(np.asarray(X_train, dtype=np.float64)),
        np.asfortranarray(np.asarray(X_test, dtype=np.float64)),
        y_train,
        y_test,
        model_features,
    )


def get_group_ids(column_names):
    """
    Given a list of the column names for X, return a list of the starting
    index of each group based on the number following the first underscore.
    The column names are expected to be in the format [cluster|kmer]_<group>_<feature>_NUM

    Note that the column names must be sorted such that all features from the same group
    are together. This is the case for the current implementation of feature generation.

    Should return an ndarry of these starting indices.
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

    # Load teh data and metadata
    data = read_feather_data(args.data)
    metadata = read_metadata(args.metadata)
    # Get the metadata columns that have more than 2 unique values
    # and more than 50 samples per category
    metadata_columns = get_metadata_columns(metadata, min_samples=args.min_samples)

    all_model_features = None
    confusion_matrix_rows = []

    # Iterate over the metadata columns
    with PdfPages(output_pdf) as pdf:
        for metadata_col in metadata_columns:
            print(f"Processing metadata column: {metadata_col}")
            print()

            X_train, X_test, y_train, y_test, model_features = merge_and_split_data(
                data,
                metadata,
                metadata_col,
                min_samples=args.min_samples,
                balanced_test=args.balanced_test,
                train_prop=args.train_prop,
            )

            # skip the column if the merge and split function returns None
            if X_train is None:
                print(
                    f"Skipping {metadata_col} as there are not enough samples after merging and filtering..."
                )
                print()
                continue

            num_classes = len(np.unique(y_train))
            print(f"Number of classes for {metadata_col}: {num_classes}")

            # set group ids based on feature names if --grouped is supplied
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

            if args.train_prop == 1:
                cm = []
                print("Not calculating test accuracy as we are not using test data...")
            else:
                # add a check to make sure there are more than 2 unique values in the predictions
                # an error can be thrown if inverse_transform gets the wrong number of columns
                # Handle cases where predictions are all of one class
                yhat = model.predict(X_test.astype(np.float64))
                if len(np.unique(yhat)) < 2:
                    print(f"Predictions for {metadata_col} are all of one class.")
                    unique_class = np.unique(yhat)[0]
                    yhat_2d = np.zeros((y_test.size, len(oh.categories_[0])))
                    yhat_2d[:, unique_class] = 1
                    y_pred = oh.inverse_transform(yhat_2d).flatten()
                    cm = confusion_matrix(y_test, y_pred)
                    accuracy = np.trace(cm) / np.sum(cm)
                    print(f"Test confusion matrix for {metadata_col}")
                    print(cm)
                    print(f"Accuracy: {accuracy:.2f}")
                else:
                    yhat_2d = np.zeros((yhat.size, len(oh.categories_[0])))
                    yhat_2d[np.arange(yhat.size), yhat] = 1
                    yhat = yhat_2d

                    try:
                        y_pred = oh.inverse_transform(yhat).flatten()
                        cm = confusion_matrix(y_test, y_pred, labels=oh.categories_[0])
                        print(f"Test confusion matrix for {metadata_col}")
                        print(cm)
                        accuracy = np.trace(cm) / np.sum(cm)
                        print(f"Accuracy: {accuracy:.2f}")
                    except Exception as e:
                        print(
                            f"Failed to transform predictions for {metadata_col}: {e}"
                        )
                        continue

            yhat_train = model.predict(X_train.astype(np.float64))
            if len(np.unique(yhat_train)) < 2:
                print(f"Train predictions for {metadata_col} are all of one class.")
                unique_class_train = np.unique(yhat_train)[0]
                yhat_train_2d = np.zeros((y_train.size, len(oh.categories_[0])))
                yhat_train_2d[:, unique_class_train] = 1
                y_train_pred = oh.inverse_transform(yhat_train_2d).flatten()
                cm_train = confusion_matrix(
                    y_train, y_train_pred, labels=oh.categories_[0]
                )
                train_accuracy = np.trace(cm_train) / np.sum(cm_train)
                print(f"Train confusion matrix for {metadata_col}")
                print(cm_train)
                print(f"Train accuracy: {train_accuracy:.2f}\n")
            else:
                yhat_train_2d = np.zeros((yhat_train.size, len(oh.categories_[0])))
                yhat_train_2d[np.arange(yhat_train.size), yhat_train] = 1
                yhat_train = yhat_train_2d

                try:
                    y_train_pred = oh.inverse_transform(yhat_train).flatten()
                    cm_train = confusion_matrix(y_train, y_train_pred)
                    print(f"Train confusion matrix for {metadata_col}")
                    print(cm_train)
                    train_accuracy = np.trace(cm_train) / np.sum(cm_train)
                    print(f"Train accuracy: {train_accuracy:.2f}\n")
                except Exception as e:
                    cm_train = []
                    print(
                        f"Failed to transform train predictions for {metadata_col}: {e}"
                    )
                    train_accuracy = 0
                    continue

            metadata_categories = oh.categories_[0]
            if args.train_prop < 1:
                append_confusion_matrix_rows(
                    confusion_matrix_rows,
                    metadata_col,
                    "test",
                    cm,
                    metadata_categories,
                )
            append_confusion_matrix_rows(
                confusion_matrix_rows,
                metadata_col,
                "train",
                cm_train,
                metadata_categories,
            )

            # extract the nonzero coefficients
            coef = model.coef_
            # get the feature names for the nonzero coefficients
            model_features = [
                f"{feature}+{category}"
                for feature in model_features
                for category in metadata_categories
            ]

            # get the names and values of the nonzero coefficients
            model_features = pd.DataFrame(model_features, columns=["feature"])
            model_features["coefficient"] = coef.toarray().flatten()
            model_features = model_features[model_features["coefficient"] != 0]

            # separate the feature names into the feature and category
            model_features["feature"] = model_features["feature"].str.split("+")
            model_features["category"] = model_features["feature"].str[1]
            model_features["feature"] = model_features["feature"].str[0]

            # gather the coefficients and categorie names into two columns grouping by feature
            model_features = (
                model_features.groupby("feature")
                .agg({"coefficient": list, "category": list})
                .reset_index()
            )
            model_features["classes"] = model_features["category"].apply(
                lambda x: "[" + ",".join(map(str, x)) + "]"
            )
            model_features["coefficients"] = model_features["coefficient"].apply(
                lambda x: "[" + ",".join(map(str, x)) + "]"
            )

            # add column for metadata category
            model_features["metadata_category"] = metadata_col

            # assuming cm can be larger than 2x2
            if args.train_prop < 1:
                if cm.shape[0] > 2:
                    accuracy = np.trace(cm) / np.sum(cm)
                    specificity = None
                    sensitivity = None
                else:
                    accuracy = (cm[0][0] + cm[1][1]) / sum(sum(row) for row in cm)
                    specificity = cm[0][0] / (cm[0][0] + cm[0][1])
                    sensitivity = cm[1][1] / (cm[1][0] + cm[1][1])
            else:
                if cm_train.shape[0] > 2:
                    accuracy = None
                    specificity = None
                    sensitivity = None
                else:
                    accuracy = None
                    specificity = cm_train[0][0] / (cm_train[0][0] + cm_train[0][1])
                    sensitivity = cm_train[1][1] / (cm_train[1][0] + cm_train[1][1])

            if args.train_prop < 1:
                if specificity is None:
                    print(f"Accuracy: {accuracy:.2f}")
                else:
                    print(
                        f"Accuracy: {accuracy:.2f}, Specificity: {specificity:.2f}, Sensitivity: {sensitivity:.2f}"
                    )
            else:
                if specificity is None:
                    print(f"Train accuracy: {train_accuracy:.2f}")
                else:
                    print(
                        f"Train Accuracy: {train_accuracy:.2f}, Train Specificity: {specificity:.2f}, Train Sensitivity: {sensitivity:.2f}"
                    )

            # add the accuracy, specificity, and sensitivity to the model features
            model_features["accuracy"] = accuracy if accuracy is not None else "NA"
            model_features["train_accuracy"] = (
                train_accuracy if train_accuracy is not None else "NA"
            )
            model_features["specificity"] = (
                specificity if specificity is not None else "NA"
            )
            model_features["sensitivity"] = (
                sensitivity if sensitivity is not None else "NA"
            )

            if args.train_prop < 1:
                out_cm = [map(str, row) for row in cm]
                out_cm = "{" + ";".join([",".join(row) for row in out_cm]) + "}"
            else:
                out_cm = None
            model_features["confusion_matrix"] = out_cm if out_cm is not None else "NA"

            model_features = model_features[
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

            # join with the larger set of model features
            if all_model_features is None:
                all_model_features = model_features
            else:
                all_model_features = pd.concat(
                    [all_model_features, model_features], axis=0
                )

            # plot the confusion matrix and save to the pdf
            # color by relative frequency
            # add numbers to the cells of the confusion matrix
            # add accuracy, specificity, and sensitivity to the title
            if args.train_prop < 1:
                plt.figure()
                plt.imshow(
                    cm / cm.sum(axis=1)[:, np.newaxis], cmap="viridis", vmin=0, vmax=1
                )
                plt.colorbar()
                for i in range(cm.shape[0]):
                    for j in range(cm.shape[1]):
                        plt.text(j, i, f"{cm[i, j]}", ha="center", va="center")
                if cm.shape[0] > 2:
                    plt.title(f"{metadata_col}\nAccuracy: {accuracy:.2f}")
                else:
                    plt.title(
                        f"{metadata_col}\nAccuracy: {accuracy:.2f}, Specificity: {specificity:.2f}, Sensitivity: {sensitivity:.2f}"
                    )
                plt.xlabel("Predicted")
                plt.ylabel("True")
                plt.xticks(range(cm.shape[1]), metadata_categories, rotation=45)
                # these need to start from the bottom
                plt.yticks(range(cm.shape[0]), metadata_categories, rotation=45)
                plt.tight_layout()
                pdf.savefig()
                plt.close()
            else:
                plt.figure()
                plt.imshow(
                    cm_train / cm_train.sum(axis=1)[:, np.newaxis],
                    cmap="viridis",
                    vmin=0,
                    vmax=1,
                )
                plt.colorbar()
                for i in range(cm_train.shape[0]):
                    for j in range(cm_train.shape[1]):
                        plt.text(j, i, f"{cm_train[i, j]}", ha="center", va="center")
                if cm_train.shape[0] > 2:
                    plt.title(f"{metadata_col}\nTrain Accuracy: {train_accuracy:.2f}")
                else:
                    plt.title(
                        f"{metadata_col}\nTrain Accuracy: {train_accuracy:.2f}, Specificity: {specificity:.2f}, Sensitivity: {sensitivity:.2f}"
                    )
                plt.xlabel("Predicted")
                plt.ylabel("True")
                plt.xticks(range(cm_train.shape[1]), metadata_categories, rotation=45)
                # these need to start from the bottom
                plt.yticks(range(cm_train.shape[0]), metadata_categories, rotation=45)
                plt.tight_layout()
                pdf.savefig()
                plt.close()

            # add a blank line
            print()

        if all_model_features is None:
            # Write a blank page to the PDF
            plt.figure()
            plt.text(
                0.5, 0.5, "No metadata columns were processed", ha="center", va="center"
            )
            plt.axis("off")
            pdf.savefig()
            plt.close()

            # Write an empty output TSV with just column names
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
            # output the nonzero coefficients to a tsv file
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
