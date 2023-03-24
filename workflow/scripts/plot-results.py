import logging

import matplotlib.pyplot as plt
import numpy as np
from astropy.convolution import convolve_fft
from astropy.io import fits
from astropy.visualization import simple_norm

from . import io

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

FIGSIZE = (5, 4)
FIGSIZE_WIDE = (13, 4)
FIGSIZE_THUMBNAIL = (1, 1)
DPI = 300
AXES_RECT = [0.11, 0.11, 0.87, 0.87]


def to_shape(a, shape):
    y_, x_ = shape
    y, x = a.shape
    y_pad = y_ - y
    x_pad = x_ - x
    return np.pad(
        a,
        ((y_pad // 2, y_pad // 2 + y_pad % 2), (x_pad // 2, x_pad // 2 + x_pad % 2)),
        mode="constant",
    )


# TODO: this plot does not seem to be all that useful, think about something better...
def plot_pixel_by_pixel_correlation(run_config, comparison_config):
    flux = run_config.flux
    flux_ref = comparison_config.flux_ground_truth

    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_axes([0.16, 0.11, 0.78, 0.8])

    non_zero = np.nonzero(flux_ref)
    lims = 0.3 * np.min(flux_ref[non_zero]), 3 * np.max(flux_ref[non_zero])

    x = np.geomspace(lims[0], lims[1], 10)
    ax.plot(x, x, color=(0.8, 0.8, 0.8))

    ax.scatter(flux_ref[non_zero], flux[non_zero], alpha=0.2)
    ax.set_xlabel("$Flux_{Ref}$")
    ax.set_ylabel("Flux")
    ax.loglog()
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_title("Flux - Pixel by Pixel Correlation")

    filename = run_config.filename_flux_correlation_image
    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_flux(run_config, comparison_config):
    flux = run_config.flux
    flux_ref = comparison_config.flux_ground_truth

    fig, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=FIGSIZE_WIDE,
        gridspec_kw={"left": 0.04, "right": 0.97, "bottom": 0.05, "top": 0.98},
    )

    norm = simple_norm(
        flux_ref,
        min_cut=0,
        max_cut=comparison_config.plot["max_cut"],
        stretch="asinh",
        asinh_a=comparison_config.plot["asinh_a"],
    )

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

    filename = run_config.filename_flux_image
    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_npred(run_config, comparison_config):
    flux = run_config.flux
    flux_ref = comparison_config.flux_ground_truth

    dataset = comparison_config.datasets_stacked
    if run_config.method == "jolideco":
        npred = (
            convolve_fft(flux + dataset["background"], dataset["psf"])
            * dataset["exposure"]
        )
    else:
        npred = convolve_fft(
            flux + (dataset["background"] * dataset["exposure"]), dataset["psf"]
        )

    npred_ref = (
        convolve_fft(flux_ref + dataset["background"], dataset["psf"])
        * dataset["exposure"]
    )

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

    filename = run_config.filename_npred_image
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


def plot_trace_lira(run_config, comparison_config):
    result = run_config.result
    result.plot_parameter_traces()
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.98)
    filename = run_config.filename_trace_image
    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_traces_jolideco(run_config, comparison_config):
    result = run_config.result
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_axes([0.16, 0.11, 0.78, 0.87])
    result.plot_trace_loss(ax=ax)

    filename = run_config.filename_trace_image
    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI // 2)


def plot_result(run_config, comparison_config):
    plot_flux(run_config=run_config, comparison_config=comparison_config)
    plot_npred(run_config=run_config, comparison_config=comparison_config)

    if run_config.method == "lira":
        plot_trace_lira(run_config=run_config, comparison_config=comparison_config)
    else:
        plot_traces_jolideco(run_config=run_config, comparison_config=comparison_config)


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
