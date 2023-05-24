import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot(input, output):
    df = pd.read_csv(input)
    df["id"] = (
        df["source_group_number"].astype(dtype=str)
        + "-"
        + df["target_group_number"].astype(dtype=str)
    )

    fig, axs = plt.subplots(nrows=2)
    fig.set_size_inches(12, 12)
    axs[0].set_title("Standard error of the mean")
    sns.barplot(data=df, x="id", y="conductivity", ci="sem", ax=axs[0])

    df = df[~df["is_neighbor"]]
    axs[1].set_title("Standard error of the mean without neighbors")
    sns.barplot(data=df, x="id", y="conductivity", ci="sem", ax=axs[1])

    plt.savefig(output)


ensemble = snakemake.input[0]  # type: ignore
output = snakemake.output[0]  # type: ignore
plot(ensemble, output)
