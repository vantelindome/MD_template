import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from dictionary import number2structure, structure2color


def plot(input, output):
    _, axs = plt.subplots(figsize=(20, 20), nrows=2, ncols=2)

    df = pd.read_csv(input)
    df = df[~df["is_neighbor"]]
    df["log_conductivity"] = np.log(np.abs(df["conductivity"]))

    for row_index in [0, 1]:
        for column_index in [0, 1]:
            ax = axs[row_index, column_index]
            ax.set_title(f"{'conductivity' if row_index==0 else 'log_conductivity'}")
            if column_index == 0:
                tmp_df = df[
                    (3 <= df["source_group_number"])
                    & (129 >= df["source_group_number"])
                    & (3 <= df["target_group_number"])
                    & (129 >= df["target_group_number"])
                ]
            else:
                tmp_df = df[
                    (130 <= df["source_group_number"])
                    & (130 <= df["target_group_number"])
                ]
            ax_heatmap = sns.heatmap(
                data=tmp_df.pivot(
                    index="source_group_number",
                    columns="target_group_number",
                    values="conductivity" if row_index == 0 else "log_conductivity",
                ),
                square=True,
                cmap="inferno_r",
                ax=ax,
            )
            for cell_index, structure in enumerate(number2structure):
                ax_heatmap.axvspan(
                    xmin=cell_index,
                    xmax=cell_index + 1,
                    ymin=0,
                    ymax=255,
                    color=structure2color[structure],
                    alpha=0.1,
                    linewidth=0,
                )
                ax_heatmap.axhspan(
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
