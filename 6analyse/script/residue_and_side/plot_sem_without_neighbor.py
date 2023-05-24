import itertools

import matplotlib.pyplot as plt
import pandas as pd


def plot(input, output):
    df = pd.read_csv(input)
    df = df[~df["is_neighbor"]]

    df.loc[df["is_neighbor"], "color"] = "red"
    df.loc[~df["is_neighbor"], "color"] = "blue"

    fig, axs = plt.subplots(ncols=3, nrows=3)
    fig.set_size_inches(24, 24)
    fig.suptitle("Standard error of mean of each pair")

    n = len(df)
    k = int(n / 9)
    for (i, j), partial_df in zip(
        itertools.product(range(3), range(3)),
        [df.loc[i : i + k - 1, :] for i in range(0, n, k)],
    ):
        axs[i][j].bar(
            partial_df.index, partial_df["conductivity"], yerr=partial_df["sem"]
        )

    plt.savefig(output)


ensemble = snakemake.input[0]  # type: ignore
output = snakemake.output[0]  # type: ignore
plot(ensemble, output)
