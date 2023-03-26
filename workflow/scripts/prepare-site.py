import shutil
import itertools
from utils import render_and_write_rst

FILENAMES_TOCTREE = {
    "Point-Sources": ["point-sources/index.rst"],
    "Disks": [f"disk{idx}/index.rst" for idx in range(1, 4)],
    "Spirals": [f"spiral{idx}/index.rst" for idx in range(1, 6)],
}


def render_index(filename):
    """Render index"""
    render_and_write_rst(
        filename=filename,
        template_name="index.rst",
        filenames_toctree=FILENAMES_TOCTREE,
    )


def render_sub_index(filename, title, sub_pages):
    """Render sub index"""
    filenames_toctree = [f"{sub_page}/index.rst" for sub_page in sub_pages]

    render_and_write_rst(
        filename=filename,
        template_name="sub-index.rst",
        filenames_toctree=filenames_toctree,
        title=title.title(),
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
            render_sub_index(
                f"results/{scenario}/{bkg_level}/index.rst",
                title=bkg_level,
                sub_pages=prefixes,
            )
            for prefix in prefixes:
                render_sub_index(
                    f"results/{scenario}/{bkg_level}/{prefix}/index.rst",
                    title=prefix,
                    sub_pages=methods,
                )

    shutil.copyfile("workflow/site/conf.py", "results/conf.py")
