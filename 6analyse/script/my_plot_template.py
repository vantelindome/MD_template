import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot(input, output):
    pass


ensemble = snakemake.input[0]  # type: ignore
output = snakemake.output[0]  # type: ignore
plot(ensemble, output)
