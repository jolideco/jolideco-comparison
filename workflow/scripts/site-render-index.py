from pathlib import Path

import yaml
from utils import render_and_write_rst

FILENAMES_TOCTREE = {
    "Point Sources": ["point1/index.rst"],
    "Point Sources & Gaussian": [f"aster{idx}/index.rst" for idx in range(1, 4)],
    "Disks": [f"disk{idx}/index.rst" for idx in range(1, 4)],
    "Spirals": [f"spiral{idx}/index.rst" for idx in range(1, 6)],
}


def read_metrics(path, methods):
    """Read metrics"""
    metrics = {}

    for method in methods:
        filename = Path(path) / f"{method}/metrics.yaml"
        with Path(filename).open("r") as f:
            metrics[method] = yaml.safe_load(f)

    return metrics


def render_index(filename):
    """Render index"""
    render_and_write_rst(
        filename=filename,
        template_name="index.rst",
        filenames_toctree=FILENAMES_TOCTREE,
    )


def render_sub_index(filename, title, sub_pages, **kwargs):
    """Render sub index"""
    filenames_toctree = [f"{sub_page}/index.rst" for sub_page in sub_pages]

    render_and_write_rst(
        filename=filename,
        template_name="sub-index.rst",
        filenames_toctree=filenames_toctree,
        title=title.title(),
        **kwargs,
    )


if __name__ == "__main__":
    render_index("results/index.rst")

    scenarios = snakemake.config["scenarios"]
    bkg_levels = snakemake.config["bkg_levels"]
    prefixes = snakemake.config["prefixes"]
    methods = snakemake.config["methods"]

    # TODO: replace by a recursive function?
    for scenario in scenarios:
        render_sub_index(
            f"results/{scenario}/index.rst",
            title=scenario,
            sub_pages=bkg_levels,
        )
        for bkg_level in bkg_levels:
            render_and_write_rst(
                filename=f"results/{scenario}/{bkg_level}/index.rst",
                template_name="sub-index-prefixes.rst",
                title=bkg_level.title(),
                methods=methods,
                prefixes=prefixes,
            )
            for prefix in prefixes:
                metrics = read_metrics(
                    path=f"results/{scenario}/{bkg_level}/{prefix}/", methods=methods
                )
                render_and_write_rst(
                    filename=f"results/{scenario}/{bkg_level}/{prefix}/index.rst",
                    template_name="sub-index-methods.rst",
                    title=prefix.title(),
                    methods=methods,
                    metrics=metrics,
                )
