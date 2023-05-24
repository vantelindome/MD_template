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

    _, ax = plt.subplots(figsize=(12, 12))
    ax.set_title("Signal pathway from FMN to J-alpha")

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

    # Calculate pathway and get edge list
    pathway = nx.shortest_path(G, source="FMN1", target="SER129", weight="length")
    edges_of_pathway = [
        set((pathway[i], pathway[i + 1])) for i in range(len(pathway) - 1)
    ]

    # Set colors of edges
    # If the edge is included in pathway, the color will be C3
    # If not, then the color will be C7
    linecollection.set_color(
        ["C3" if set((u, v)) in edges_of_pathway else "C7" for (u, v) in G.edges()]
    )

    plt.axis("off")
    plt.savefig(output)


ensemble = snakemake.input[0]  # type: ignore
output = snakemake.output[0]  # type: ignore
plot(ensemble, output)
