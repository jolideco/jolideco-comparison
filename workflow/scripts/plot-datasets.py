import logging

import matplotlib.pyplot as plt
from astropy.visualization import simple_norm
from plot import AXES_RECT, DPI, FIGSIZE
from utils import to_shape

from .plot import AXES_RECT, DPI, FIGSIZE, FIGSIZE_THUMBNAIL, FIGSIZE_WIDE

log = logging.getLogger(__name__)


def plot_exposure(comparison_config):
    """Plot exposure image."""
    exposure = comparison_config.datasets_stacked["exposure"]

    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_axes(AXES_RECT)

    im = ax.imshow(exposure, origin="lower")

    ax.set_xlabel("x / pix")
    ax.set_ylabel("y / pix")
    fig.colorbar(im, ax=ax)

    filename = comparison_config.filename_exposure_image
    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_psf(comparison_config):
    """Plot PSF image."""
    psf = comparison_config.datasets_stacked["psf"]

    shape = comparison_config.datasets_stacked["counts"].shape
    psf = to_shape(psf, shape=shape)

    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_axes(AXES_RECT)

    norm = simple_norm(
        psf,
        min_cut=0,
        stretch="asinh",
        asinh_a=0.1,
    )

    im = ax.imshow(psf, norm=norm, origin="lower")

    ax.set_xlabel("x / pix")
    ax.set_ylabel("y / pix")
    fig.colorbar(im, ax=ax)

    filename = comparison_config.filename_psf_image
    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_counts(comparison_config):
    """Plot counts image."""
    counts = comparison_config.datasets_stacked["counts"]

    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_axes(AXES_RECT)

    norm = simple_norm(
        counts,
        min_cut=0,
        stretch="asinh",
        asinh_a=0.1,
    )

    im = ax.imshow(counts, norm=norm, origin="lower", interpolation="None")

    ax.set_xlabel("x / pix")
    ax.set_ylabel("y / pix")
    fig.colorbar(im, ax=ax)

    filename = comparison_config.filename_counts_image
    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_background(comparison_config):
    """Plot background image."""
    background = comparison_config.datasets_stacked["background"]

    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_axes(AXES_RECT)

    im = ax.imshow(background, origin="lower")

    ax.set_xlabel("x / pix")
    ax.set_ylabel("y / pix")
    fig.colorbar(im, ax=ax)

    filename = comparison_config.filename_background_image
    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_dataset(comparison_config):
    """Plot dataset."""
    plot_psf(comparison_config=comparison_config)
    plot_counts(comparison_config=comparison_config)
    plot_exposure(comparison_config=comparison_config)
    plot_background(comparison_config=comparison_config)
