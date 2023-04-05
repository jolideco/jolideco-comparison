from pathlib import Path

import numpy as np
import yaml
from utils import read_sub_config, render_and_write_rst


def get_deconvolver_configuration(config):
    """Get configuration"""
    from jolideco.core import MAPDeconvolver
    from pylira.core import LIRADeconvolver

    if config["method"] == "jolideco":
        deconvolver = MAPDeconvolver(**config["deconvolver"])
    elif config["method"] == "lira":
        deconvolver = LIRADeconvolver(**config["deconvolver"])
    else:
        raise ValueError(f"Unknown method: {config['method']}")

    return str(deconvolver) + "\n"


def get_model_configuration(config):
    """Get model configuration"""
    from jolideco.models import FluxComponents

    if config["method"] == "jolideco":
        config["components"]["flux"]["flux_upsampled"] = np.zeros((128, 128))
        components = FluxComponents.from_dict(config["components"])
        return str(components) + "\n"

    elif config["method"] == "lira":
        return
    else:
        raise ValueError(f"Unknown method: {config['method']}")


def render_summary(filename, filename_result, config):
    """Render summary"""
    configuration = get_deconvolver_configuration(config=config)
    model_configuration = get_model_configuration(config=config)

    title = config["name"].replace("-", " ").title()

    with (filename.parent / "metrics.yaml").open("r") as f:
        metrics = yaml.safe_load(f)

    render_and_write_rst(
        filename=filename,
        template_name="summary.rst",
        title=title,
        configuration=configuration,
        model_configuration=model_configuration,
        metrics=metrics,
        filename_result=filename_result,
    )


if __name__ == "__main__":
    config = read_sub_config(snakemake.input[0], method=snakemake.wildcards.method)

    filename_reult = Path(snakemake.output[4]).name
    render_summary(
        Path(snakemake.output[0]), config=config, filename_result=filename_reult
    )
