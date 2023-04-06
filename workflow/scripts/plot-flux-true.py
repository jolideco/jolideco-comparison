import numpy as np
from astropy.io import fits
from plot import plot_flux_thumbnail
from utils import read_config

if __name__ == "__main__":
    flux_ref = fits.getdata(snakemake.input[0]).astype(np.float32)

    config = read_config(snakemake.input.filename_config)

    filename = snakemake.output.filename_image_flux_ref
    plot_flux_thumbnail(flux_ref, config=config, filename=filename)
