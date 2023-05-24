import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from dictionary import number2chain


def plot(input, output):
    df = pd.read_csv(input)
    df["log_conductivity"] = np.log(np.abs(df["conductivity"]))
    for side in ["source", "target"]:
        df[f"{side}_group_chain"] = df["source_group_number"].map(
            lambda n: number2chain[n - 1]
        )
    # Remove flavin and inter monomer
    # Leave only intra monomer
    df = df[df["source_group_chain"] == df["target_group_chain"]]

    chainA = df[df["source_group_chain"] == "A"]
    chainB = df[df["target_group_chain"] == "B"]

    for side in ["source", "target"]:
        chainA[f"{side}_group_number"] = chainA[f"{side}_group_number"] - 2
        chainB[f"{side}_group_number"] = chainB[f"{side}_group_number"] - 129

        chainA[f"{side}_group_id"] = (
            chainA[f"{side}_group_number"].astype(str) + chainA[f"{side}_group_name"]
        )
        chainB[f"{side}_group_id"] = (
            chainB[f"{side}_group_number"].astype(str) + chainB[f"{side}_group_name"]
        )

    merge_df = pd.merge(
        left=chainA,
        right=chainB,
        on=["source_group_id", "target_group_id"],
        suffixes=["_A", "_B"],
    )

    plt.figure(figsize=(10, 10))
    conductivity = "log_conductivity"
    sns.jointplot(
        data=merge_df,
        x=f"{conductivity}_A",
        y=f"{conductivity}_B",
        hue="is_neighbor_B",
    )
    plt.savefig(output)


ensemble = snakemake.input[0]  # type: ignore
output = snakemake.output[0]  # type: ignore
plot(ensemble, output)
