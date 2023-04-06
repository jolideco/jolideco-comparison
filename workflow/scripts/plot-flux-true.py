import numpy as np
from astropy.io import fits
from plot import plot_flux_thumbnail
from utils import read_config

if __name__ == "__main__":
    flux_ref = fits.getdata(snakemake.input[0]).astype(np.float32)

    path = f"config/{snakemake.wildcards.scenario}/bg1/chandra.yaml"
    config = read_config(path)

    plot_flux_thumbnail(flux_ref, config=config, filename=snakemake.output[0])
