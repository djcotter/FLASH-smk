from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import OneHotEncoder

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import adelie as ad

import numpy as np

import pyarrow.feather as feather
import pandas as pd
from math import floor

import argparse

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
        help=(
            "Proportion of the data to use for training. Grabs this proportion "
            "from the smallest class and then evenly samples that number from "
            "all other classes."
        ),
    )
    return parser.parse_args()


def read_feather_data(file_path):
    return feather.read_feather(file_path)


def read_metadata(file_path):
    metadata = pd.read_table(file_path)
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


def balance_and_split_data_into_folds(
    data, metadata, metadata_col, min_samples=50, n_folds=10
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
        return None

    merged_data = merged_data[merged_data[metadata_col].isin(classes_to_keep)]

    # Balance the data by keeping the same number of samples for each class
    num_to_keep = class_counts.min()
    if pd.isna(num_to_keep):
        return None

    balanced_indices = (
        merged_data.groupby(metadata_col)
        .apply(lambda x: x.sample(n=num_to_keep, replace=False).index)
        .explode()
    )
    balanced_data = merged_data.loc[balanced_indices]

    # Shuffle the balanced data to ensure randomness
    balanced_data = balanced_data.sample(frac=1, random_state=42).reset_index(drop=True)

    # Split the balanced data into n_folds ensuring unique test data for each fold
    folds = []
    fold_size = floor(len(balanced_data) / n_folds)

    for fold_idx in range(n_folds):
        start_idx = fold_idx * fold_size
        end_idx = (
            start_idx + fold_size if fold_idx < n_folds - 1 else len(balanced_data)
        )

        test_indices = balanced_data.iloc[start_idx:end_idx].index
        train_indices = balanced_data.index.difference(test_indices)

        X_train = balanced_data.drop(["sample_name", metadata_col], axis=1).loc[
            train_indices
        ]
        X_test = balanced_data.drop(["sample_name", metadata_col], axis=1).loc[
            test_indices
        ]
        y_train = balanced_data[metadata_col].loc[train_indices].to_numpy()
        y_test = balanced_data[metadata_col].loc[test_indices].to_numpy()
        model_features = X_train.columns

        folds.append(
            (
                np.asfortranarray(X_train),
                np.asfortranarray(X_test),
                y_train,
                y_test,
                model_features,
            )
        )

    return folds


def train_adelie_model(X_train, y_train, n_threads=1):
    oh = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
    y_train2 = oh.fit_transform(y_train[:, np.newaxis])

    model = ad.GroupElasticNet(solver="cv_grpnet", family="multinomial")
    model.fit(
        X_train.astype(np.float64), y_train2.astype(np.float64), n_threads=n_threads
    )

    return model, oh


def main():
    args = parse_args()
    output_prefix = args.output_prefix
    output_pdf = output_prefix + "_confusion_matrices.pdf"
    output_coef = output_prefix + "_nonzero_coefficients.tsv"

    # Load teh data and metadata
    data = read_feather_data(args.data)
    metadata = read_metadata(args.metadata)
    # Get the metadata columns that have more than 2 unique values
    # and more than 50 samples per category
    metadata_columns = get_metadata_columns(metadata, min_samples=args.min_samples)

    all_model_features = None

    # Iterate over the metadata columns
    with PdfPages(output_pdf) as pdf:
        for metadata_col in metadata_columns:
            print(f"Processing metadata column: {metadata_col}")
            print()

            folds = balance_and_split_data_into_folds(
                data,
                metadata,
                metadata_col,
                min_samples=args.min_samples,
                n_folds=10,
            )

            if folds is None:
                print(
                    f"Skipping {metadata_col} as there are not enough samples after merging and filtering..."
                )
                print()
                continue

            for fold_idx, (
                X_train,
                X_test,
                y_train,
                y_test,
                model_features,
            ) in enumerate(folds):
                print(
                    f"Processing fold {fold_idx + 1} for metadata column: {metadata_col}"
                )
                print()

                try:
                    model, oh = train_adelie_model(
                        X_train, y_train, n_threads=args.n_threads
                    )
                except Exception as e:
                    print(
                        f"Failed to train model for {metadata_col}, fold {fold_idx + 1}: {e}"
                    )
                    continue

                # Handle cases where predictions are all of one class
                yhat = model.predict(X_test.astype(np.float64))
                if len(np.unique(yhat)) < 2:
                    print(f"Predictions for {metadata_col} are all of one class.")
                    unique_class = np.unique(yhat)[0]
                    yhat_2d = np.zeros((y_test.size, yhat.max() + 1))
                    yhat_2d[:, unique_class] = 1
                    y_pred = oh.inverse_transform(yhat_2d).flatten()
                    cm = confusion_matrix(y_test, y_pred)
                    accuracy = np.trace(cm) / np.sum(cm)
                    print(f"Test confusion matrix for {metadata_col}")
                    print(cm)
                    print(f"Accuracy: {accuracy:.2f}")
                else:
                    yhat_2d = np.zeros((yhat.size, yhat.max() + 1))
                    yhat_2d[np.arange(yhat.size), yhat] = 1
                    yhat = yhat_2d

                    try:
                        y_pred = oh.inverse_transform(yhat).flatten()
                        cm = confusion_matrix(y_test, y_pred)
                        print(f"Test confusion matrix for {metadata_col}")
                        print(cm)
                    except Exception as e:
                        print(
                            f"Failed to transform predictions for {metadata_col}: {e}"
                        )
                        continue

                yhat_train = model.predict(X_train.astype(np.float64))
                if len(np.unique(yhat_train)) < 2:
                    print(f"Train predictions for {metadata_col} are all of one class.")
                    unique_class_train = np.unique(yhat_train)[0]
                    yhat_train_2d = np.zeros((y_train.size, yhat_train.max() + 1))
                    yhat_train_2d[:, unique_class_train] = 1
                    y_train_pred = oh.inverse_transform(yhat_train_2d).flatten()
                    cm_train = confusion_matrix(y_train, y_train_pred)
                    train_accuracy = np.trace(cm_train) / np.sum(cm_train)
                    print(
                        f"Train confusion matrix for {metadata_col}, fold {fold_idx + 1}"
                    )
                    print(cm_train)
                    print(f"Train accuracy: {train_accuracy:.2f}\n")
                else:
                    yhat_train_2d = np.zeros((yhat_train.size, yhat_train.max() + 1))
                    yhat_train_2d[np.arange(yhat_train.size), yhat_train] = 1
                    yhat_train = yhat_train_2d

                    try:
                        y_train_pred = oh.inverse_transform(yhat_train).flatten()
                        cm_train = confusion_matrix(y_train, y_train_pred)
                        print(
                            f"Train confusion matrix for {metadata_col}, fold {fold_idx + 1}"
                        )
                        print(cm_train)
                        train_accuracy = np.trace(cm_train) / np.sum(cm_train)
                        print(f"Train accuracy: {train_accuracy:.2f}\n")
                    except Exception as e:
                        cm_train = []
                        print(
                            f"Failed to transform train predictions for {metadata_col}, fold {fold_idx + 1}: {e}"
                        )
                        train_accuracy = 0
                        continue

                # Extract nonzero coefficients
                coef = model.coef_
                metadata_categories = oh.categories_[0]
                model_features = [
                    f"{feature}+{category}"
                    for feature in model_features
                    for category in metadata_categories
                ]

                model_features = pd.DataFrame(model_features, columns=["feature"])
                model_features["coefficient"] = coef.toarray().flatten()
                model_features = model_features[model_features["coefficient"] != 0]

                model_features["feature"] = model_features["feature"].str.split("+")
                model_features["category"] = model_features["feature"].str[1]
                model_features["feature"] = model_features["feature"].str[0]

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

                model_features["metadata_category"] = metadata_col
                model_features["fold_index"] = fold_idx + 1

                accuracy = (
                    np.trace(cm) / np.sum(cm)
                    if cm.shape[0] > 2
                    else (cm[0][0] + cm[1][1]) / np.sum(cm)
                )
                specificity = (
                    cm[0][0] / (cm[0][0] + cm[0][1]) if cm.shape[0] == 2 else None
                )
                sensitivity = (
                    cm[1][1] / (cm[1][0] + cm[1][1]) if cm.shape[0] == 2 else None
                )

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

                out_cm = [map(str, row) for row in cm]
                out_cm = "{" + ";".join([",".join(row) for row in out_cm]) + "}"
                model_features["confusion_matrix"] = (
                    out_cm if out_cm is not None else "NA"
                )

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
                        "fold_index",
                    ]
                ]

                if all_model_features is None:
                    all_model_features = model_features
                else:
                    all_model_features = pd.concat(
                        [all_model_features, model_features], axis=0
                    )

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
                        f"{metadata_col}, Fold {fold_idx + 1}\nAccuracy: {accuracy:.2f}"
                    )
                else:
                    plt.title(
                        f"{metadata_col}, Fold {fold_idx + 1}\nAccuracy: {accuracy:.2f}, Specificity: {specificity:.2f}, Sensitivity: {sensitivity:.2f}"
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


if __name__ == "__main__":
    main()
