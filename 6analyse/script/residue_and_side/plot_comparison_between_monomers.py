import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from dictionary import number2chain


def add_confidence_interval(ax, mean, std):
    for i in range(1, 4):
        color = "C" + str(i + 1)
        ax.axhspan(
            ymin=mean + (i - 1) * std,
            ymax=mean + i * std,
            facecolor=color,
            alpha=0.2,
            label=f"mu + {str(i)}sigma",
        )
        ax.axhspan(
            ymin=mean - (i - 1) * std,
            ymax=mean - i * std,
            facecolor=color,
            alpha=0.2,
        )
        ax.legend()


def plot(input, output):
    df = pd.read_csv(input)
    for side in ["source", "target"]:
        df[f"{side}_group_chain"] = df[f"{side}_group_number"].map(
            lambda n: number2chain[n - 1]
        )
        df = df[df[f"{side}_group_chain"] != "FMN"]

    # Get only interactions amoung each monomer
    df = df[df["source_group_chain"] == df["target_group_chain"]]

    chains = set(number2chain)
    chains.remove("FMN")
    chains = list(chains)
    # chainA, chainB = [df[df[f"source_group_chain"] == chain] for chain in chains]
    chainA = df[df["source_group_chain"] == "A"]
    chainB = df[df["source_group_chain"] == "B"]
    chainA["pair_id"] = (
        chainA["source_group_number"].astype(str)
        + "-"
        + chainA["target_group_number"].astype(str)
    )
    tmp_source_id = chainB["source_group_number"] - 127
    tmp_target_id = chainB["target_group_number"] - 127
    chainB["pair_id"] = tmp_source_id.astype(str) + "-" + tmp_target_id.astype(str)

    compared_df = pd.merge(
        left=chainA,
        right=chainB,
        how="inner",
        on=["pair_id"],
        suffixes=["_A", "_B"],
    )
    compared_df["conductivity"] = (
        compared_df["conductivity_A"] - compared_df["conductivity_B"]
    )
    compared_df["log_conductivity"] = np.log(np.abs(compared_df["conductivity"]))

    fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(12, 12))
    fig.suptitle("Difference of conductivities between monomers")

    ax_heat = sns.heatmap(
        data=compared_df.pivot(
            index="source_group_number_A",
            columns="target_group_number_A",
            values="log_conductivity",
        ),
        square=True,
        cmap="inferno_r",
        ax=axs[0][0],
    )
    ax_heat.set_title("Log of abs of difference of heat conductivity between monomers")

    ax_scatter = sns.scatterplot(
        x=compared_df.index,
        y=compared_df.conductivity,
        hue=compared_df.is_neighbor_A,
        marker="+",
        ax=axs[1][0],
    )
    add_confidence_interval(
        ax=ax_scatter,
        mean=compared_df["conductivity"].mean(),
        std=compared_df["conductivity"].std(),
    )
    ax_scatter.set_title("Diff of conductivities between monomers")
    ax_scatter.xaxis.set_visible(False)

    ax_hist = sns.histplot(
        data=compared_df,
        y="log_conductivity",
        ax=axs[1][1],
        kde=True,
        hue="is_neighbor_A",
    )
    ax_hist.axhline(
        y=np.log(compared_df["conductivity"].mean()),
        color="C3",
        label="average of difference",
    )
    ax_hist.legend()
    ax_hist.set_title("Population of log of abs of conductivity")

    compared_df = compared_df[~compared_df["is_neighbor_A"]]
    axs[0][1].set_title("Diff of conductivities between monomers\nwithout neighbors")
    axs[0][1].scatter(
        x=compared_df["pair_id"], y=compared_df["conductivity"], marker="+"
    )
    add_confidence_interval(
        ax=axs[0][1],
        mean=compared_df["conductivity"].mean(),
        std=compared_df["conductivity"].std(),
    )
    axs[0][1].xaxis.set_visible(False)

    plt.tight_layout()
    plt.savefig(output)


ensemble = snakemake.input[0]  # type: ignore
output = snakemake.output[0]  # type: ignore
plot(ensemble, output)
