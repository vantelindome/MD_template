import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from dictionary import number2structure, structure2color


def plot(input, output):
    _, axs = plt.subplots(figsize=(30, 15), ncols=2)

    df = pd.read_csv(input)
    df = df[~df["is_neighbor"]]
    df["log_conductivity"] = np.log(np.abs(df["conductivity"]))

    for ax, conductivity in zip(axs, ["conductivity", "log_conductivity"]):
        ax = sns.heatmap(
            data=df.pivot(
                index="source_group_number",
                columns="target_group_number",
                values=conductivity,
            ),
            cmap="inferno_r",
            square=True,
            ax=ax,
        )
        ax.set_title(conductivity)

        for cell_index, structure in enumerate(number2structure):
            ax.axvspan(
                xmin=cell_index,
                xmax=cell_index + 1,
                ymin=0,
                ymax=255,
                color=structure2color[structure],
                alpha=0.1,
                linewidth=0,
            )
            ax.axhspan(
                ymin=cell_index,
                ymax=cell_index + 1,
                xmin=0,
                xmax=255,
                color=structure2color[structure],
                alpha=0.1,
                linewidth=0,
            )

    plt.savefig(output)


ensemble = snakemake.input[0]  # type: ignore
output = snakemake.output[0]  # type: ignore
plot(ensemble, output)
