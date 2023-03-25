from pathlib import Path

import docutils.core
from jinja2 import Environment, FileSystemLoader
from utils import read_sub_config


def rst_to_html(rst):
    """Convert reStructuredText to HTML"""
    return docutils.core.publish_parts(source=rst, writer_name="html")["html_body"]


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


def render_rst_jinja(config):
    """Render HTML"""
    environment = Environment(loader=FileSystemLoader("workflow/report/"))
    template = environment.get_template("summary-run-template.rst")

    configuration = get_deconvolver_configuration(config=config)
    title = config["name"].replace("-", " ").title()
    return template.render(title=title, configuration=configuration)


if __name__ == "__main__":
    config = read_sub_config(snakemake.input[0], method=snakemake.wildcards.method)

    rst_rendered = render_rst_jinja(config=config)
    html_rendered = rst_to_html(rst_rendered)

    with Path(snakemake.output[0]).open("w") as f:
        f.write(html_rendered)
