import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from dictionary import chain2color, number2chain

number2chain_include_fmn = number2chain
number2chain_include_fmn[0] = "A"
number2chain_include_fmn[1] = "B"


def plot(input, output):
    df = pd.read_csv(input)

    for side in ["source", "target"]:
        df[f"{side}_group_name"] = df[f"{side}_group_number"].astype(str)
        df[f"{side}_group_chain"] = df[f"{side}_group_number"].map(
            lambda n: number2chain_include_fmn[n - 1]
        )

    df = df[df["source_group_chain"] == df["target_group_chain"]]
    chainA = df[df["source_group_chain"] == "A"]
    chainB = df[df["source_group_chain"] == "B"]
    chainA["pair_id"] = (
        chainA["source_group_number"].astype(str)
        + "-"
        + chainA["target_group_number"].astype(str)
    )

    # FMN: 1, 2
    # chain A: 3-129
    # chain B: 130-255
    # Now we treat FMN 1 as chain A and treat FMN 2 as chain B
    # A: 1, 3-129
    # B: 2, 130-255
    tmp_source_id = chainB["source_group_number"] - 127
    tmp_target_id = chainB["target_group_number"] - 127
    tmp_source_id = tmp_source_id.where(tmp_source_id > 0, 1)
    tmp_target_id = tmp_target_id.where(tmp_target_id > 0, 1)

    chainB["pair_id"] = tmp_source_id.astype(str) + "-" + tmp_target_id.astype(str)

    compared_df = pd.merge(
        left=chainA,
        right=chainB,
        on=["pair_id"],
        suffixes=["_A", "_B"],
    )
    compared_df["diff_conductivity"] = (
        compared_df["conductivity_A"] - compared_df["conductivity_B"]
    )
    compared_df["mean_conductivity"] = (
        compared_df["conductivity_A"] + compared_df["conductivity_B"] / 2
    )
    compared_df["length"] = 1 / compared_df["mean_conductivity"]
    compared_df["log_diff_conductivity"] = np.log(
        np.abs(compared_df["diff_conductivity"])
    )

    avg = compared_df["diff_conductivity"].mean()
    std = compared_df["diff_conductivity"].std()
    compared_df["in_1_sigma"] = np.abs(compared_df["diff_conductivity"] - avg) < std
    compared_df["in_2_sigma"] = np.abs(compared_df["diff_conductivity"] - avg) < 2 * std
    compared_df["in_3_sigma"] = np.abs(compared_df["diff_conductivity"] - avg) < 3 * std
    for side in ["source", "target"]:
        compared_df[f"{side}_group_name_A"] = compared_df[
            f"{side}_group_name_A"
        ].astype(str)

    G = nx.from_pandas_edgelist(
        df=compared_df,
        source="source_group_name_A",
        target="target_group_name_A",
        edge_attr=True,
    )

    nodes = (
        pd.DataFrame(
            {
                "node": [str(i) for i in range(1, 129)],
                "color": [
                    chain2color[number2chain_include_fmn[i - 1]]
                    if i != 1
                    else chain2color["FMN"]
                    for i in range(1, 129)
                ],
            }
        )
        .set_index("node")
        .to_dict(orient="index")
    )

    _, ax = plt.subplots(figsize=(12, 12))
    pos = nx.spring_layout(G, seed=1)

    nx.set_node_attributes(G, nodes)
    nx.draw_networkx_nodes(
        G, pos, node_color=[d["color"] for (_, d) in G.nodes(data=True)], ax=ax
    )

    nx.draw_networkx_labels(G, pos, ax=ax)

    width = lambda x: np.abs(x) * 1e2  # np.exp(x)

    linecollection = nx.draw_networkx_edges(
        G,
        pos,
        width=[width(d["diff_conductivity"]) for (_, _, d) in G.edges(data=True)],
        ax=ax,
    )

    # Controller of edge colors
    # Change me if you would like to change edge color
    interval2color = {
        "in_1_sigma": "C2",
        "in_2_sigma": "C3",
        "in_3_sigma": "C4",
        "outside": "C5",
    }

    linecollection.set_colors(
        [
            interval2color["in_1_sigma"]
            if d["in_1_sigma"]
            else interval2color["in_2_sigma"]
            if d["in_2_sigma"]
            else interval2color["in_3_sigma"]
            if d["in_3_sigma"]
            else interval2color["outside"]
            for (_, _, d) in G.edges(data=True)
        ]
    )
    opacity = 0.4
    linecollection.set_alpha(opacity)
    ax.legend(
        handles=[
            mpatches.Patch(color=v, label=k, alpha=opacity)
            for k, v in interval2color.items()
        ]
    )
    ax.set_title("Comparison between monomers")
    plt.axis("off")
    plt.savefig(output)


ensemble = snakemake.input[0]  # type: ignore
output = snakemake.output[0]  # type: ignore
plot(ensemble, output)
