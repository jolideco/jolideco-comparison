import logging

import matplotlib.pyplot as plt
from astropy.visualization import simple_norm
from plot import AXES_RECT, DPI, FIGSIZE
from utils import read_datasets_all, stack_datasets, to_shape

log = logging.getLogger(__name__)


def plot_exposure(exposure, filename):
    """Plot exposure image."""
    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_axes(AXES_RECT)

    im = ax.imshow(exposure, origin="lower")

    ax.set_xlabel("x / pix")
    ax.set_ylabel("y / pix")
    fig.colorbar(im, ax=ax)

    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_psf(psf, filename):
    """Plot PSF image."""
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

    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_counts(counts, filename):
    """Plot counts image."""
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

    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_background(background, filename):
    """Plot background image."""
    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_axes(AXES_RECT)

    im = ax.imshow(background, origin="lower")

    ax.set_xlabel("x / pix")
    ax.set_ylabel("y / pix")
    fig.colorbar(im, ax=ax)

    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


if __name__ == "__main__":
    datasets = read_datasets_all(
        scenario=snakemake.wildcards.scenario,
        bkg_level=snakemake.wildcards.bkg_level,
        prefix=snakemake.wildcards.prefix,
    )

    stacked = stack_datasets(datasets)

    plot_counts(counts=stacked["counts"], filename=snakemake.output[0])
    plot_exposure(exposure=stacked["exposure"], filename=snakemake.output[1])
    plot_background(background=stacked["background"], filename=snakemake.output[2])

    stacked["psf"] = to_shape(stacked["psf"], shape=stacked["counts"].shape)
    plot_psf(psf=stacked["psf"], filename=snakemake.output[3])
