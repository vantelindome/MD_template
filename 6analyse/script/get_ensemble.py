import re

import numpy as np
import pandas as pd


def get_ensemble(input_files, output):
    dfs = [
        pd.read_table(
            file,
            delim_whitespace=True,
            header=None,
            names=["source_group_id", "target_group_id", "conductivity"],
            dtype={
                "source_group_id": str,
                "target_group_id": str,
                "conductivity": float,
            },
        )
        for file in input_files
    ]

    groupby = pd.concat(dfs).groupby(
        by=["source_group_id", "target_group_id"], sort=True, as_index=False
    )

    ensemble_df = groupby.mean()
    for operator in ["std", "var", "count", "sem"]:
        method = getattr(groupby, operator)
        ensemble_df[operator] = method()["conductivity"]

    # Add columns source_group_number, source_group_name,
    # target_group_number, and target_group_name
    for side in ["source", "target"]:
        ensemble_df[f"{side}_group_number"] = ensemble_df[f"{side}_group_id"].map(
            lambda col: int(re.sub("^0+", "", col.split("_")[0]))
        )
        ensemble_df[f"{side}_group_name"] = ensemble_df[f"{side}_group_id"].map(
            lambda col: col.split("_")[1]
            if "-" not in col.split("_")[1]
            else col.split("_")[1].split("-")[0]
        )
    # Add columns if the pair next to each other
    ensemble_df["is_neighbor"] = (
        np.abs(ensemble_df["source_group_number"] - ensemble_df["target_group_number"])
        == 1
    )
    # Add columns to use plot network
    ensemble_df["length"] = 1 / ensemble_df["conductivity"].abs()

    ensemble_df.to_csv(output, index=False)


files = snakemake.input  # type: ignore
output = snakemake.output[0]  # type: ignore


get_ensemble(files, output)
