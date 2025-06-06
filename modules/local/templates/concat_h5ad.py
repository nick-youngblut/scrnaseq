#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import scanpy as sc, anndata as ad, pandas as pd
from pathlib import Path
import platform


def read_samplesheet(samplesheet):
    df = pd.read_csv(samplesheet)
    df.set_index("sample")

    # samplesheet may contain replicates, when it has,
    # group information from replicates and collapse with commas
    # only keep unique values using set()
    df = df.groupby(["sample"]).agg(lambda column: ",".join(set(column.astype(str))))

    return df

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.
    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.
    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str

def dump_versions():
    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "scanpy": sc.__version__,
        }
    }

    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions))

if __name__ == "__main__":

    # Open samplesheet as dataframe
    df_samplesheet = read_samplesheet("${samplesheet}")

    # find all h5ad and append to list, keeping sample names
    adatas = []
    sample_names = []
    for path in Path(".").rglob("*.h5ad"):
        adatas.append(sc.read_h5ad(path))
        sample_names.append(path.stem.replace("_matrix", ""))

    # concatenate using outer join on features
    adata = ad.concat(adatas, merge="outer", label="sample", keys=sample_names)

    # merge with data.frame, on sample information
    adata.obs = adata.obs.join(df_samplesheet, on="sample", how="left").astype(str)
    adata.write_h5ad("${meta.id}.h5ad")

    print("Wrote h5ad file to {}".format("${meta.id}.h5ad"))

    # dump versions
    dump_versions()
