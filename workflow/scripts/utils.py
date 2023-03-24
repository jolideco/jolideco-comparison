import copy
import logging
from pathlib import Path

import numpy as np
import yaml
from astropy.io import fits

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

INSTRUMENTS = {
    "chandra": ["chandra"],
    "xmm": ["xmm"],
    "joint": ["chandra", "xmm"],
}

DATA_PATH = {
    "chandra": "sims_chandra_gauss_fwhm4710_128x128",
    "xmm": "sims_xmm_gauss_fwhm14130_128x128",
}

FILE_PATTERN_CHANDRA = {
    "counts": "chandra_*_sim*_{bkg_level}_{name}*iter*.fits",
    "psf": "chandra_gauss_fwhm4710_128x128_psf_33x33.fits",
    "flux": "chandra_gauss_fwhm4710_128x128_flux_33x33.fits",
    "npred": "chandra_gauss_fwhm4710_128x128_npred_33x33.fits",
}

FILE_PATTERN_XMM = {
    "counts": "xmm_*_sim*_{bkg_level}_{name}*iter*.fits",
    "psf": "xmm_gauss_fwhm14130_128x128_psf_63x63.fits",
    "flux": "xmm_gauss_fwhm14130_128x128_psf_63x63.fits",
    "npred": "chandra_gauss_fwhm4710_128x128_npred_33x33.fits",
}

FILE_PATTERN = {
    "chandra": FILE_PATTERN_CHANDRA,
    "xmm": FILE_PATTERN_XMM,
}


def to_shape(data, shape):
    """Pad array to shape"""
    y_pad = shape[0] - data.shape[0]
    x_pad = shape[1] - data.shape[1]

    pad_width = (
        (y_pad // 2, y_pad // 2 + y_pad % 2),
        (x_pad // 2, x_pad // 2 + x_pad % 2),
    )
    return np.pad(data, pad_width, mode="constant")


def read_config(filename):
    """Read config"""
    log.info(f"Reading {filename}")

    with Path(filename).open("r") as f:
        config = yaml.safe_load(f)

    return config


def get_instrument_and_idx(filename):
    """Get instrument and idx from filename"""
    parts = filename.stem.split("_")
    instrument, idx = parts[0], parts[-1]
    return instrument, idx


def get_filenames_counts(instrument, bkg_level, name):
    """Find files"""
    path = Path(f"data") / DATA_PATH[instrument]

    if name == "point-sources":
        name = ""

    pattern = FILE_PATTERN[instrument]["counts"].format(bkg_level=bkg_level, name=name)

    return list(path.glob(pattern))


def read_dataset(filename_counts, filename_psf):
    """Read single dataset"""
    counts = fits.getdata(filename_counts).astype(np.float32)
    psf = fits.getdata(filename_psf).astype(np.float32)

    exposure = np.ones_like(counts)
    background = np.zeros_like(counts)

    return {
        "counts": counts,
        "psf": psf,
        "exposure": exposure,
        "background": background,
    }


def read_datasets(filenames_counts):
    """Read multiple datasets"""
    datasets = {}

    for filename in filenames_counts:
        instrument, idx = get_instrument_and_idx(filename)
        name = f"{instrument}-{idx}"

        filename_psf = filename.parent / FILE_PATTERN[instrument]["psf"]

        dataset = read_dataset(filename_counts=filename, filename_psf=filename_psf)
        datasets[name] = dataset

    return datasets


def read_datasets_all(prefix, bkg_level, name):
    """Read all datasets"""
    datasets_all = {}

    for instrument in INSTRUMENTS[prefix]:
        filename_counts = get_filenames_counts(
            instrument=instrument, bkg_level=bkg_level, name=name
        )
        datasets = read_datasets(filenames_counts=filename_counts)
        datasets_all.update(datasets)

    return datasets_all


def stack_datasets(datasets):
    """Stack datasets"""
    datasets = copy.deepcopy(datasets)
    stacked = datasets.popitem()[1]

    psf_shapes = np.array([dataset["psf"].shape for dataset in datasets.values()])
    psf_shape = tuple(psf_shapes.max(axis=0))

    for dataset in datasets.values():
        for key, value in dataset.items():
            if key == "psf":
                value = to_shape(value, psf_shape)
            stacked[key] += value

    return stacked
