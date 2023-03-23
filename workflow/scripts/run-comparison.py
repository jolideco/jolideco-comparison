import logging
from pathlib import Path

import numpy as np
import yaml
from astropy.io import fits

# from jolideco.core import MAPDeconvolver
# from jolideco.models import FluxComponents

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

DATA_PATH = {
    "chandra": "sims_chandra_gauss_fwhm4710_128x128",
    "xmm": "sims_xmm_gauss_fwhm14130_128x128",
    "joint": "sims_chandra_gauss_fwhm4710_128x128,sims_xmm_gauss_fwhm14130_128x128",
}

PSF_PATH = {
    "chandra": "chandra_gauss_fwhm4710_128x128_psf_33x33.fits",
    "xmm": "xmm_gauss_fwhm14130_128x128_psf_63x63.fits",
}


def read_config(filename):
    """Read config"""
    log.info(f"Reading {filename}")

    with Path(filename).open("r") as f:
        config = yaml.safe_load(f)

    return config


def read_dataset(filename_counts):
    """Read single dataset"""
    counts = fits.getdata(filename_counts)

    instrument = str(filename_counts.stem.split("_")[0])
    idx = str(filename_counts.stem.split("_")[-1])

    psf = fits.getdata(filename_counts.parent / PSF_PATH[instrument])

    exposure = np.ones_like(counts)
    background = np.zeros_like(counts)

    name = f"{instrument}-{idx}"

    data = {
        "counts": counts,
        "psf": psf,
        "exposure": exposure,
        "background": background,
    }
    return {name: data}


def read_datasets(config):
    """Read datasets"""
    datasets = {}

    filename_counts = []

    for sub_path in DATA_PATH[config["prefix"]].split(","):
        path = Path(f"data/{sub_path}")

        if config["name"] != "point-sources":
            pattern = f"*_{config['bkg_level']}_{config['name']}_*.fits"
        else:
            pattern = f"*_sim01_{config['bkg_level']}_*.fits"

        filename_counts.extend(path.glob(pattern))

    print(filename_counts)

    for filename in filename_counts:
        dataset = read_dataset(filename)
        datasets.update(dataset)

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

    print(datasets)
    for run in config["runs"]:
        if run["method"] == "jolideco":
            run_jolideco(datasets)
        elif run["method"] == "lira":
            run_pylira(datasets)
        else:
            raise ValueError(f"Unknown run: {run['method']}")


if __name__ == "__main__":
    config = read_config(snakemake.input[0])

    run_comparison(config)
