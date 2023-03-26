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


def render_summary(filename, config):
    """Render summary"""
    configuration = get_deconvolver_configuration(config=config)
    title = config["name"].replace("-", " ").title()

    render_and_write_rst(
        filename=filename,
        template_name="summary.rst",
        title=title,
        configuration=configuration,
    )


if __name__ == "__main__":
    config = read_sub_config(snakemake.input[0], method=snakemake.wildcards.method)

    render_summary(snakemake.output[0], config=config)
