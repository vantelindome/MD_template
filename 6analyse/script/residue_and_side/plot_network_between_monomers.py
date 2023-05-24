import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from dictionary import chain2color, number2chain


def is_inter_monomer(row):
    a = [1, *range(3, 130)]
    b = [2, *range(130, 256)]
    u = row["source_group_number"]
    v = row["target_group_number"]
    if u in a:
        if v in b:
            return True
        else:
            return False
    elif v in b:
        if u in a:
            return True
        else:
            return False


def plot(input, output):
    df = pd.read_csv(input)

    add_chain_name = (
        lambda n: f"{n}{number2chain[n-1]}"
        if n <= 128
        else f"{n-128}{number2chain[n-1]}"
    )
    for side in ["source", "target"]:
        df[f"{side}_group_name"] = df[f"{side}_group_number"].map(add_chain_name)

    df["is_inter_monomer"] = [is_inter_monomer(r) for _, r in df.iterrows()]

    G = nx.from_pandas_edgelist(
        df=df, source="source_group_name", target="target_group_name", edge_attr=True
    )

    nodes = (
        pd.DataFrame(
            {
                "node": [add_chain_name(i) for i in range(1, 256)],
                "color": [chain2color[number2chain[i - 1]] for i in range(1, 256)],
            }
        )
        .set_index("node")
        .to_dict(orient="index")
    )

    _, ax = plt.subplots(figsize=(12, 12))
    ax.set_title("Network between monomers")
    pos = nx.spring_layout(G, seed=1)

    nx.set_node_attributes(G, nodes)
    nx.draw_networkx_nodes(
        G, pos, node_color=[d["color"] for (_, d) in G.nodes(data=True)], ax=ax
    )

    nx.draw_networkx_labels(G, pos, ax=ax)

    linecollection = nx.draw_networkx_edges(
        G,
        pos,
        width=[
            d["conductivity"] * 1000 / df["conductivity"].sum()
            for (u, v, d) in G.edges(data=True)
        ],
        ax=ax,
    )
    linecollection.set_colors(
        ["C3" if d["is_inter_monomer"] else "C7" for (_, _, d) in G.edges(data=True)]
    )

    plt.axis("off")
    plt.savefig(output)


ensemble = snakemake.input[0]  # type: ignore
output = snakemake.output[0]  # type: ignore
plot(ensemble, output)
