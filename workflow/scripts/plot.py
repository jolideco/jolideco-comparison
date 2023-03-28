import logging

import matplotlib.pyplot as plt
from astropy.visualization import simple_norm

FIGSIZE = (5, 4)
FIGSIZE_WIDE = (13, 4)
FIGSIZE_THUMBNAIL = (1, 1)

DPI = 300
AXES_RECT = [0.11, 0.11, 0.87, 0.87]

log = logging.getLogger(__name__)


def plot_flux_thumbnail(flux, filename, config):
    """Plot flux thumbnail"""
    fig = plt.figure(figsize=FIGSIZE_THUMBNAIL)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()

    kwargs = config["plot"]["flux"]["norm"]
    kwargs["max_cut"] = kwargs["max_cut"]
    norm = simple_norm(flux, **kwargs)

    ax.imshow(flux, norm=norm, origin="lower")
    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=DPI)
