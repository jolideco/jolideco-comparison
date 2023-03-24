import logging

import matplotlib.pyplot as plt
import numpy as np
from astropy.convolution import convolve_fft
from astropy.visualization import simple_norm
from utils import (
    read_config,
    read_datasets_all,
    read_deconvolution_result,
    read_flux_ref,
    stack_datasets,
    to_shape,
)

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

FIGSIZE = (5, 4)
FIGSIZE_WIDE = (13, 4)
FIGSIZE_THUMBNAIL = (1, 1)
DPI = 300
AXES_RECT = [0.11, 0.11, 0.87, 0.87]


def plot_flux(flux, flux_ref, filename, config):
    """Plot flux and residuals"""
    fig, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=FIGSIZE_WIDE,
        gridspec_kw={"left": 0.04, "right": 0.97, "bottom": 0.05, "top": 0.98},
    )

    norm = simple_norm(flux_ref, **config["plot"]["flux"]["norm"])

    im = axes[0].imshow(flux, norm=norm, origin="lower")

    axes[0].set_xlabel("x / pix")
    axes[0].set_ylabel("y / pix")
    axes[0].set_title("Flux")
    fig.colorbar(im, ax=axes[0])

    im = axes[1].imshow(flux_ref, norm=norm, origin="lower")

    axes[1].set_xlabel("x / pix")
    axes[1].set_ylabel("y / pix")
    axes[1].set_title("$Flux_{Ref}$")
    fig.colorbar(im, ax=axes[1])

    residual = (flux - flux_ref) / np.sqrt(flux_ref)

    im = axes[2].imshow(residual, origin="lower", vmin=-2, vmax=2, cmap="RdBu")

    axes[2].set_xlabel("x / pix")
    axes[2].set_ylabel("y / pix")
    axes[2].set_title("Residuals")
    cbar = fig.colorbar(im, ax=axes[2])
    cbar.ax.set_ylabel("(Flux - $Flux_{Ref}$) / $\sqrt{Flux_{Ref}}$")

    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_npred(npred, npred_ref, filename):
    """Plot npred and residuals"""
    fig, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=FIGSIZE_WIDE,
        gridspec_kw={"left": 0.04, "right": 0.97, "bottom": 0.05, "top": 0.98},
    )

    norm = simple_norm(
        npred_ref,
        min_cut=0,
        stretch="asinh",
        asinh_a=0.1,
    )

    im = axes[0].imshow(npred, norm=norm, origin="lower")

    axes[0].set_xlabel("x / pix")
    axes[0].set_ylabel("y / pix")
    axes[0].set_title("$NPred$")
    fig.colorbar(im, ax=axes[0])

    im = axes[1].imshow(npred_ref, norm=norm, origin="lower")

    axes[1].set_xlabel("x / pix")
    axes[1].set_ylabel("y / pix")
    axes[1].set_title("$Npred_{Ref}$")
    fig.colorbar(im, ax=axes[1])

    residual = (npred - npred_ref) / np.sqrt(npred_ref)

    im = axes[2].imshow(residual, origin="lower", vmin=-2, vmax=2, cmap="RdBu")

    axes[2].set_xlabel("x / pix")
    axes[2].set_ylabel("y / pix")
    axes[2].set_title("Residuals")
    cbar = fig.colorbar(im, ax=axes[2])
    cbar.ax.set_ylabel("(NPred - $NPred_{Ref}$) / $\\sqrt{NPred_{Ref}}$")

    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_flux_thumbnail(run_config, comparison_config):
    flux = run_config.flux

    fig = plt.figure(figsize=FIGSIZE_THUMBNAIL)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()

    flux_ref = comparison_config.flux_ground_truth

    norm = simple_norm(
        flux_ref,
        min_cut=0,
        max_cut=comparison_config.plot["max_cut"],
        stretch="asinh",
        asinh_a=comparison_config.plot["asinh_a"],
    )

    ax.imshow(flux, norm=norm, origin="lower")
    filename = run_config.filename_flux_thumbnail
    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_flux_reference_thumbnail(comparison_config):
    flux_ref = comparison_config.flux_ground_truth

    fig = plt.figure(figsize=FIGSIZE_THUMBNAIL)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()

    norm = simple_norm(
        flux_ref,
        min_cut=0,
        max_cut=comparison_config.plot["max_cut"],
        stretch="asinh",
        asinh_a=comparison_config.plot["asinh_a"],
    )

    ax.imshow(flux_ref, norm=norm, origin="lower")
    filename = comparison_config.filename_ground_truth_thumbnail
    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_counts_thumbnail(comparison_config):
    counts = comparison_config.datasets_stacked["counts"]

    fig = plt.figure(figsize=FIGSIZE_THUMBNAIL)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()

    norm = simple_norm(
        counts,
        min_cut=0,
        max_cut=5 * comparison_config.plot["max_cut"],
        stretch="asinh",
        asinh_a=0.1,
    )

    ax.imshow(counts, norm=norm, origin="lower")
    filename = comparison_config.filename_counts_thumbnail
    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_trace_lira(result, filename):
    """Plot trace of LIRA parameters."""
    result.plot_parameter_traces()
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.98)
    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_traces_jolideco(result, filename):
    """Plot trace of Jolideco parameters."""
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_axes([0.16, 0.11, 0.78, 0.87])
    result.plot_trace_loss(ax=ax)

    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI // 2)


def plot_exposure(comparison_config):
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
    plot_psf(comparison_config=comparison_config)
    plot_counts(comparison_config=comparison_config)
    plot_exposure(comparison_config=comparison_config)
    plot_background(comparison_config=comparison_config)


if __name__ == "__main__":
    kwargs = {
        "name": snakemake.wildcards.name,
        "bkg_level": snakemake.wildcards.bkg_level,
        "prefix": snakemake.wildcards.prefix,
        "method": snakemake.wildcards.method,
    }

    result = read_deconvolution_result(**kwargs)
    datasets = read_datasets_all(**kwargs)
    dataset = stack_datasets(datasets=datasets)

    if "pylira" in kwargs["method"]:
        plot_trace_lira(result=result, filename=snakemake.output[2])
        flux = result.flux

        npred = convolve_fft(
            flux + (dataset["background"] * dataset["exposure"]), dataset["psf"]
        )

    if "jolideco" in kwargs["method"]:
        plot_traces_jolideco(result=result, filename=snakemake.output[2])
        flux = result.flux_numpy
        npred = (
            convolve_fft(flux + dataset["background"], dataset["psf"])
            * dataset["exposure"]
        )

    path = "config/{name}/{bkg_level}/{prefix}.yaml".format(**kwargs)
    config = read_config(path)

    flux_ref = read_flux_ref(**kwargs)
    npred_ref = read_npred_ref(**kwargs)

    plot_flux(
        flux=flux,
        flux_ref=flux_ref,
        filename=snakemake.output[0],
        config=config,
    )

    plot_npred(
        npred=npred,
        npred_ref=npred_ref,
        filename=snakemake.output[1],
        config=config,
    )
