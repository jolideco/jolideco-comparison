import copy

import numpy as np
from utils import read_config, read_datasets, stack_datasets

RANDOM_STATE = np.random.RandomState(7362)

N_ITER_DEBUG = 5


def prepare_datasets_lira(datasets):
    """Prepare datasets for LIRA"""
    stacked = stack_datasets(datasets=datasets)
    # TODO: this is a possible bug in Pylira...
    stacked["background"] = stacked["background"] / stacked["exposure"]
    return stacked


def prepare_datasets_jolideco(datasets):
    """Prepare datasets for jolideco"""
    datasets = copy.deepcopy(datasets)

    for dataset in datasets.values():
        dataset["psf"] = {"flux": dataset["psf"]}

    return datasets


def get_flux_init(datasets, oversample=10.0):
    """Get flux init"""
    stacked = stack_datasets(datasets=datasets)

    flux = stacked["counts"] / stacked["exposure"] - stacked["background"]

    flux_mean = np.nanmean(np.clip(flux, 0, np.inf))

    flux_init = RANDOM_STATE.gamma(oversample * flux_mean, size=flux.shape) / oversample
    return flux_init.astype(np.float32)


def run_jolideco(datasets, config, debug):
    """Run jolideco"""
    from jolideco.core import MAPDeconvolver
    from jolideco.models import FluxComponents, NPredCalibration, NPredCalibrations

    flux_init = get_flux_init(datasets=datasets)

    if config["components"]["flux"].get("upsampling_factor", 1) > 1:
        flux_init = flux_init.repeat(2, axis=0).repeat(2, axis=1)

    config["components"]["flux"]["flux_upsampled"] = flux_init
    components = FluxComponents.from_dict(config["components"])

    deconvolver = MAPDeconvolver(**config["deconvolver"])

    if debug:
        deconvolver.n_epochs = N_ITER_DEBUG

    datasets = prepare_datasets_jolideco(datasets=datasets)

    calibrations = NPredCalibrations()

    for name in datasets:
        calibration = NPredCalibration(background_norm=1)
        calibration.shift_xy.requires_grad = False
        calibrations[name] = calibration

    result = deconvolver.run(
        datasets=datasets, components=components, calibrations=calibrations
    )
    return result


def run_pylira(datasets, config, debug):
    """Run LIRA"""
    from pylira import LIRADeconvolver

    dataset = prepare_datasets_lira(datasets=datasets)

    dataset["flux_init"] = get_flux_init(datasets=datasets)

    deconvolver = LIRADeconvolver(**config["deconvolver"])

    if debug:
        deconvolver.n_iter_max = N_ITER_DEBUG

    result = deconvolver.run(data=dataset)
    return result.reduce_to_mean_std()


def run_deconvolution(datasets, config_run, debug):
    """Run deconvolution"""
    if config_run["method"] == "jolideco":
        return run_jolideco(datasets, config=config_run, debug=debug)
    elif config_run["method"] == "lira":
        return run_pylira(datasets, config=config_run, debug=debug)
    else:
        raise ValueError(f"Unknown method: {config_run['method']}")


if __name__ == "__main__":
    debug = snakemake.config["debug"]
    config = read_config(snakemake.input.config)

    datasets = read_datasets(
        filenames_counts=snakemake.input.filenames_counts,
        filenames_psf=snakemake.input.filenames_psf,
    )

    for idx, config_run in enumerate(config["runs"]):
        result = run_deconvolution(
            datasets=datasets, config_run=config_run, debug=debug
        )
        result.write(snakemake.output[idx])
