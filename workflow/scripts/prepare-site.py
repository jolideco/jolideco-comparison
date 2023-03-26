import shutil
import itertools
from utils import render_and_write_rst


def render_index(filename, config):
    """Render index"""
    filenames_toctree = []

    for parts in itertools.product(
        config["scenarios"], config["bkg_levels"], config["prefixes"], config["methods"]
    ):
        scenario, bkg_level, prefix, method = parts
        filenames_toctree.append(
            f"{scenario}/{bkg_level}/{prefix}/{method}/summary.rst"
        )

    render_and_write_rst(
        filename=filename,
        template_name="index.rst",
        filenames_toctree=filenames_toctree,
    )


if __name__ == "__main__":
    render_index("results/index.rst", config=snakemake.config)

    shutil.copyfile("workflow/site/conf.py", "results/conf.py")
