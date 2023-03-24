import copy

import numpy as np
from utils import read_config, read_datasets_all, stack_datasets

DEBUG = True

RANDOM_STATE = np.random.RandomState(7362)


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


def get_flux_init(datasets):
    """Get flux init"""
    stacked = stack_datasets(datasets=datasets)

    flux = stacked["counts"] / stacked["exposure"] - stacked["background"]

    flux_init = np.clip(flux, 0, np.inf)

    # flux_init = RANDOM_STATE.gamma(flux_mean, size=flux.shape).astype(np.float32)
    return flux_init.astype(np.float32)


def run_jolideco(datasets, config):
    """Run jolideco"""
    from jolideco.core import MAPDeconvolver
    from jolideco.models import FluxComponents

    datasets = prepare_datasets_jolideco(datasets=datasets)

    config["components"]["flux"]["flux_init"] = get_flux_init(datasets=datasets)
    components = FluxComponents.from_dict(config["components"])

    deconvolver = MAPDeconvolver(**config["deconvolver"])
    result = deconvolver.run(datasets=datasets, components=components)
    return result


def run_pylira(datasets, config):
    """Run LIRA"""
    from pylira import Deconvolver

    dataset = prepare_datasets_lira(datasets=datasets)

    dataset["flux_init"] = get_flux_init(datasets=datasets)

    deconvolver = Deconvolver(**config["deconvolver"])
    result = deconvolver.run(data=dataset)
    return result


def run_deconvolution(datasets, config_run):
    """Run deconvolution"""
    if config_run["method"] == "jolideco":
        return run_jolideco(datasets, config=config_run)
    elif config_run["method"] == "lira":
        return run_pylira(datasets, config=config_run)
    else:
        raise ValueError(f"Unknown method: {config_run['method']}")


if __name__ == "__main__":
    config = read_config(snakemake.input[0])

    datasets = read_datasets_all(
        prefix=snakemake.wildcards.prefix,
        bkg_level=snakemake.wildcards.bkg_level,
        name=snakemake.wildcards.name,
    )

    for config_run in config["runs"]:
        if config_run["name"] == snakemake.wildcards.method:
            break

    result = run_deconvolution(datasets=datasets, config_run=config_run)

    result.write(snakemake.output[0])
