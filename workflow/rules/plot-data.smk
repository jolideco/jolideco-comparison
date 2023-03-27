rule plot_data:
    input:
        "results/{scenario}/{bkg_level}/{prefix}/{method}/{scenario}-{bkg_level}-{method}-{prefix}-result.fits.gz"
    output:
        "results/{scenario}/images/flux-true-thumbnail.png",
    log:
        "logs/plot-data/{scenario}-{bkg_level}-{method}-{prefix}.log"
    script:
        "../scripts/plot-results.py"