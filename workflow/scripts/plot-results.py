import logging

import matplotlib.pyplot as plt
import numpy as np
import yaml
from astropy.convolution import convolve_fft
from astropy.visualization import simple_norm
from plot import DPI, FIGSIZE_THUMBNAIL, FIGSIZE_WIDE, plot_flux_thumbnail
from skimage import metrics
from utils import (
    read_config,
    read_datasets_all,
    read_deconvolution_result,
    read_flux_ref,
    stack_datasets,
)

log = logging.getLogger(__name__)

METRICS = {
    "MSE": metrics.mean_squared_error,
    "SSI": metrics.structural_similarity,
}


def compute_metrics(flux, flux_ref):
    """Compute metrics"""
    results = {}

    for name, metric in METRICS.items():
        if name == "SSI":
            data_range = flux_ref.max() - flux_ref.min()
            value = metric(flux, flux_ref, data_range=data_range)
        else:
            value = metric(flux, flux_ref)

        results[name] = float(value)

    return results


def plot_flux(flux, flux_ref, filename, config):
    """Plot flux and residuals"""
    fig, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=FIGSIZE_WIDE,
        gridspec_kw={"left": 0.04, "right": 0.97, "bottom": 0.05, "top": 0.98},
    )

    kwargs = config["plot"]["flux"]["norm"]
    kwargs["max_cut"] = 5 * kwargs["max_cut"]
    norm = simple_norm(flux_ref, **kwargs)

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

    with np.errstate(invalid="ignore", divide="ignore"):
        residual = np.nan_to_num((flux - flux_ref) / np.sqrt(flux_ref))

    im = axes[2].imshow(residual, origin="lower", vmin=-2, vmax=2, cmap="RdBu")

    axes[2].set_xlabel("x / pix")
    axes[2].set_ylabel("y / pix")
    axes[2].set_title("Residuals")
    cbar = fig.colorbar(im, ax=axes[2])
    cbar.ax.set_ylabel("(Flux - $Flux_{Ref}$) / $\sqrt{Flux_{Ref}}$")

    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_npred(npred, npred_ref, filename, config):
    """Plot npred and residuals"""
    fig, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=FIGSIZE_WIDE,
        gridspec_kw={"left": 0.04, "right": 0.97, "bottom": 0.05, "top": 0.98},
    )

    norm = simple_norm(flux_ref, **config["plot"]["npred"]["norm"])

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

    with np.errstate(invalid="ignore", divide="ignore"):
        residual = np.nan_to_num((npred - npred_ref) / np.sqrt(npred_ref))

    im = axes[2].imshow(residual, origin="lower", vmin=-2, vmax=2, cmap="RdBu")

    axes[2].set_xlabel("x / pix")
    axes[2].set_ylabel("y / pix")
    axes[2].set_title("Residuals")
    cbar = fig.colorbar(im, ax=axes[2])
    cbar.ax.set_ylabel("(NPred - $NPred_{Ref}$) / $\\sqrt{NPred_{Ref}}$")

    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)


def plot_counts_thumbnail(counts, filename, config):
    """Plot counts thumbnail"""
    fig = plt.figure(figsize=FIGSIZE_THUMBNAIL)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    kwargs = config["plot"]["flux"]["norm"]
    kwargs["max_cut"] = 5 * kwargs["max_cut"]
    norm = simple_norm(flux_ref, **kwargs)

    ax.imshow(counts, norm=norm, origin="lower")
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


def npred_pylira(flux, dataset):
    """Npred for Pylira"""
    npred = convolve_fft(
        flux + (dataset["background"] * dataset["exposure"]), dataset["psf"]
    )
    return npred


def npred_jolideco(flux, dataset):
    """Npred for Jolideco"""
    npred = (
        convolve_fft(flux + dataset["background"], dataset["psf"]) * dataset["exposure"]
    )
    return npred


if __name__ == "__main__":
    kwargs = {
        "scenario": snakemake.wildcards.scenario,
        "bkg_level": snakemake.wildcards.bkg_level,
        "prefix": snakemake.wildcards.prefix,
    }

    result = read_deconvolution_result(snakemake.input[0])
    datasets = read_datasets_all(**kwargs)
    dataset = stack_datasets(datasets=datasets)

    flux_ref = read_flux_ref(scenario=kwargs["scenario"])

    if "pylira" in snakemake.wildcards.method:
        plot_trace_lira(result=result, filename=snakemake.output[2])
        flux = result.posterior_mean_from_trace
        npred = npred_pylira(flux=flux, dataset=dataset)
        npred_ref = npred_pylira(flux=flux_ref, dataset=dataset)

    if "jolideco" in snakemake.wildcards.method:
        plot_traces_jolideco(result=result, filename=snakemake.output[2])
        flux = result.flux_total
        npred = npred_jolideco(flux=flux, dataset=dataset)
        npred_ref = npred_jolideco(flux=flux_ref, dataset=dataset)

    metrics = compute_metrics(flux=flux, flux_ref=flux_ref)

    log.info(f"Writing {snakemake.output[4]}")

    with open(snakemake.output[4], "w") as f:
        yaml.dump(metrics, f)

    path = "config/{scenario}/{bkg_level}/{prefix}.yaml".format(**kwargs)
    config = read_config(path)

    plot_flux(
        flux=flux,
        flux_ref=flux_ref,
        filename=snakemake.output[0],
        config=config,
    )

    plot_flux_thumbnail(
        flux=flux,
        filename=snakemake.output[3],
        config=config,
    )

    plot_npred(
        npred=npred,
        npred_ref=npred_ref,
        filename=snakemake.output[1],
        config=config,
    )
