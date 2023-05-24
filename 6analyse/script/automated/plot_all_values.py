import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot(input, output):
    df = pd.read_csv(input)
    df = df[df["count"] != 1]
    _, axd = plt.subplot_mosaic([["conductivity", "var"], ["std", "sem"]])

    for key, ax in axd.items():
        sns.histplot(
            data=df,
            x=key,
            hue="is_neighbor",
            ax=ax,
            log_scale=True,
        )
    plt.tight_layout()
    plt.savefig(output)


ensemble = snakemake.input[0]  # type: ignore
output = snakemake.output[0]  # type: ignore
plot(ensemble, output)
