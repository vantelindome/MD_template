import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def plot(input, output):
    df = pd.read_csv(input)
    df["conductivity"] = np.log(np.abs(df["conductivity"]))

    fig, ax = plt.subplots()
    fig.set_size_inches(12, 12)
    ax = sns.heatmap(
        data=df.pivot(
            index="source_group_number",
            columns="target_group_number",
            values="conductivity",
        ),
        square=True,
        cmap="inferno_r",
        ax=ax,
    )
    ax.set_title("Pairwised log of conductivity")
    plt.savefig(output)


ensemble = snakemake.input[0]  # type: ignore
output = snakemake.output[0]  # type: ignore
plot(ensemble, output)
