import shutil
import itertools
from utils import render_and_write_rst

FILENAMES_TOCTREE = {
    "Point-Sources Scenario": ["point-sources/index.rst"],
    "Disks Scenario": [f"disk{idx}/index.rst" for idx in range(1, 4)],
    "Spirals Scenario": [f"spiral{idx}/index.rst" for idx in range(1, 6)],
}


def render_index(filename):
    """Render index"""
    render_and_write_rst(
        filename=filename,
        template_name="index.rst",
        filenames_toctree=FILENAMES_TOCTREE,
    )


def render_sub_index(filename, config, scenario):
    """Render index"""
    filenames_toctree = []

    for parts in itertools.product(
        config["bkg_levels"], config["prefixes"], config["methods"]
    ):
        bkg_level, prefix, method = parts
        filenames_toctree.append(f"{bkg_level}/{prefix}/{method}/summary.rst")

    render_and_write_rst(
        filename=filename,
        template_name="sub-index.rst",
        filenames_toctree=filenames_toctree,
        title=f"{scenario.title()}",
    )


if __name__ == "__main__":
    render_index(
        "results/index.rst",
    )

    for scenario in snakemake.config["scenarios"]:
        render_sub_index(
            f"results/{scenario}/index.rst",
            config=snakemake.config,
            scenario=scenario,
        )

    shutil.copyfile("workflow/site/conf.py", "results/conf.py")
