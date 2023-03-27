from plot import plot_flux_thumbnail
from utils import read_config, read_flux_ref

if __name__ == "__main__":
    flux_ref = read_flux_ref(scenario=snakemake.wildcards.scenario)

    path = f"config/{snakemake.wildcards.scenario}/bg1/chandra.yaml"
    config = read_config(path)

    plot_flux_thumbnail(flux_ref, config=config, filename=snakemake.output[0])
