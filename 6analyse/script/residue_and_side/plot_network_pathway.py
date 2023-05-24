import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from dictionary import chain2color, number2chain, number2structure


def plot(input, output):
    df = pd.read_csv(input)

    add_chain_name = (
        lambda n: f"{n}{number2chain[n-1]}"
        if n <= 2
        else f"{n-2}{number2chain[n-1]}"
        if n <= 129
        else f"{n-128}{number2chain[n-1]}"
    )
    for side in ["source", "target"]:
        df[f"{side}_group_name"] = df[f"{side}_group_number"].map(add_chain_name)

    G = nx.from_pandas_edgelist(
        df=df, source="source_group_name", target="target_group_name", edge_attr=True
    )
    pathway = nx.shortest_path(G, source="1FMN", target="126A", weight="length")

    nodes = (
        pd.DataFrame(
            {
                "node": [add_chain_name(i) for i in range(1, 256)],
                "color": [chain2color[number2chain[i - 1]] for i in range(1, 256)],
                "secondary_structure": [number2structure[i - 1] for i in range(1, 256)],
            }
        )
        .set_index("node")
        .to_dict(orient="index")
    )

    _, ax = plt.subplots(figsize=(12, 12))
    ax.set_title("Signal pathway from FMN to J-alpha")
    pos = nx.spring_layout(G, seed=1)

    nx.set_node_attributes(G, nodes)
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
    linecollection.set_color(
        ["C3" if u in pathway and v in pathway else "C7" for (u, v) in G.edges()]
    )

    ax.legend()

    plt.axis("off")
    plt.savefig(output)


ensemble = snakemake.input[0]  # type: ignore
output = snakemake.output[0]  # type: ignore
plot(ensemble, output)
