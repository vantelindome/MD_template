import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot(input, output):
    df = pd.read_csv(input)
    df = df[df["count"] != 1]
    _ = plt.subplot()
    sns.pairplot(data=df, vars=["conductivity", "var", "std", "sem"], hue="is_neighbor")
    plt.savefig(output)


ensemble = snakemake.input[0]  # type: ignore
output = snakemake.output[0]  # type: ignore
plot(ensemble, output)
