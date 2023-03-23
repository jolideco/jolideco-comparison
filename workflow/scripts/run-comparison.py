import logging

import numpy as np
import yaml
from astropy.io import fits

# from jolideco.core import MAPDeconvolver
# from jolideco.models import FluxComponents

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def read_config(filename):
    """Read config"""
    log.info(f"Reading {filename}")

    with filename.open("r") as f:
        config = yaml.safe_load(f)

    return config


def read_dataset():
    """Read single dataset"""
    counts = fits.getdata()
    psf = fits.getdata()

    exposure = np.ones_like(counts)
    background = np.zeros_like(counts)

    return {
        "counts": counts,
        "psf": psf,
        "exposure": exposure,
        "background": background,
    }


def read_datasets():
    """Read datasets"""
    datasets = {}

    for name in []:
        datasets[name] = read_dataset()

    return datasets


def run_jolideco(datasets):
    components = FluxComponents.from_dict()

    deco = MAPDeconvolver()
    result = deco.run(datasets=datasets, components=components)
    result.write()


def run_pylira():
    pass


def run_comparison(config):
    datasets = read_datasets(config)
    run_jolideco(datasets)
    run_pylira(datasets)


if __name__ == "__main__":
    print("Hello")
    config = read_config(snakemake.input[0])
    print(config)

    # run_comparison(config)
