import math

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from dictionary import chain2color, number2chain


def number2residue(df, n):
    tmp_df = df[df["source_group_number"] == n]
    if len(tmp_df) != 0:
        return tmp_df["source_group_name"].unique()[0]
    else:
        return df[df["target_group_number"] == n]["target_group_name"].unique()[0]


def my_color(min, max, n_region, value):
    width = (max - min) / n_region
    i = 1
    if math.isnan(value):
        return "white"
    while True:
        if value <= min + width * i:
            color = f"C{i}"
            break
        else:
            i += 1
    return color


def my_label2color(min, max, n_region):
    width = (max - min) / n_region
    return {
        f"{min + width * i} =< value =< {min + width * (i + 1)}": f"C{i+1}"
        for i in range(n_region)
    }


def plot(input, output):
    df = pd.read_csv(input)

    for side in ["source", "target"]:
        df[f"{side}_group_name_with_number"] = df[f"{side}_group_name"] + df[
            f"{side}_group_number"
        ].astype(str)

    G = nx.from_pandas_edgelist(
        df=df,
        source="source_group_name_with_number",
        target="target_group_name_with_number",
        edge_attr=True,
    )

    nodes = (
        pd.DataFrame(
            {
                "node": [number2residue(df, i) + str(i) for i in range(1, 256)],
                "color": [chain2color[number2chain[i - 1]] for i in range(1, 256)],
            }
        )
        .set_index("node")
        .to_dict(orient="index")
    )
    nx.set_node_attributes(G, nodes)

    pos = nx.spring_layout(G, seed=1)

    _, axs = plt.subplots(figsize=(24, 12), ncols=2)
    for ax, quantity in zip(axs, ["std", "var"]):
        ax.set_title(f"{quantity} of conductivity")
        ax.axis("off")

        nx.draw_networkx_nodes(
            G, pos, node_color=[d["color"] for (_, d) in G.nodes(data=True)], ax=ax
        )

        nx.draw_networkx_labels(G, pos, ax=ax)

        ave = df["conductivity"].sum()
        linecollection = nx.draw_networkx_edges(
            G,
            pos,
            width=[d["conductivity"] * 1000 / ave for (u, v, d) in G.edges(data=True)],
            ax=ax,
        )

        quantity_min = df[quantity].min()
        quantity_max = df[quantity].max()
        n_region = 4
        opacity = 0.6

        linecollection.set_color(
            [
                my_color(quantity_min, quantity_max, n_region, d["var"])
                for (_, _, d) in G.edges(data=True)
            ]
        )
        linecollection.set_alpha(opacity)

        ax.legend(
            handles=[
                mpatches.Patch(color=v, label=k, alpha=opacity)
                for k, v in my_label2color(quantity_min, quantity_max, n_region).items()
            ]
        )
    plt.tight_layout()
    plt.savefig(output)


ensemble = snakemake.input[0]  # type: ignore
output = snakemake.output[0]  # type: ignore
plot(ensemble, output)
