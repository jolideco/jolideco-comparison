import copy
import logging
from pathlib import Path

import numpy as np
import yaml
from astropy.io import fits
from jinja2 import Environment, FileSystemLoader

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def get_bkg_level(filename):
    """Get background level"""
    filename = Path(filename).stem

    if "bg1" in filename:
        return 0.01
    elif "bg2" in filename:
        return 0.1
    elif "bg3" in filename:
        return 1.0

    raise ValueError(f"Cannot get backgrounnd level from: {filename}")


def to_path_list(filenames):
    """Convert filenames to path list"""
    return [Path(filename) for filename in filenames]


def to_shape(data, shape):
    """Pad array to shape"""
    y_pad = shape[0] - data.shape[0]
    x_pad = shape[1] - data.shape[1]

    pad_width = (
        (y_pad // 2, y_pad // 2 + y_pad % 2),
        (x_pad // 2, x_pad // 2 + x_pad % 2),
    )
    return np.pad(data, pad_width, mode="constant")


def read_deconvolution_result(filename):
    """Read deconvolution result"""
    from jolideco import MAPDeconvolverResult
    from pylira import LIRADeconvolverResult

    if "jolideco" in str(filename):
        return MAPDeconvolverResult.read(filename)
    else:
        return LIRADeconvolverResult.read(filename)


def read_config(filename):
    """Read config"""
    log.info(f"Reading {filename}")

    with Path(filename).open("r") as f:
        config = yaml.safe_load(f)

    return config


def read_sub_config(filename, method):
    """Read sub config"""
    config = read_config(filename)

    for config_run in config["runs"]:
        if config_run["name"] == method:
            return config_run
    else:
        raise ValueError(f"Method {method} not found in {filename}")


def get_instrument_and_idx(filename):
    """Get instrument and idx from filename"""
    parts = Path(filename).stem.split("_")
    instrument, idx = parts[0], parts[-1]
    return instrument, idx


def read_npred_ref(instrument, name):
    """Read reference npred"""
    filename = get_filenames(
        instrument=instrument, bkg_level="", name=name, quantity="npred"
    )[0]

    flux_ref = fits.getdata(filename).astype(np.float32)

    return flux_ref


def read_dataset(filename_counts, filename_psf):
    """Read single dataset"""
    counts = fits.getdata(filename_counts).astype(np.float32)
    psf = fits.getdata(filename_psf).astype(np.float32)

    exposure = np.ones_like(counts)

    bkg_level = get_bkg_level(filename_counts)
    background = bkg_level * np.ones_like(counts)

    return {
        "counts": counts,
        "psf": psf,
        "exposure": exposure,
        "background": background,
    }


def read_datasets(filenames_counts, filenames_psf):
    """Read multiple datasets"""
    datasets = {}

    for filename_counts, filename_psf in zip(filenames_counts, filenames_psf):
        instrument, idx = get_instrument_and_idx(filename_counts)
        name = f"{instrument}-{idx}"

        dataset = read_dataset(
            filename_counts=filename_counts, filename_psf=filename_psf
        )
        datasets[name] = dataset

    return datasets


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


def render_and_write_rst(filename, template_name, **kwargs):
    """Render RST"""
    environment = Environment(loader=FileSystemLoader("workflow/site/templates"))
    environment.globals.update(zip=zip)

    template = environment.get_template(template_name)

    rst_rendered = template.render(**kwargs)

    with Path(filename).open("w") as f:
        f.write(rst_rendered)
