rule plot_results:
    input:
        "results/{scenario}/{bkg_level}/{prefix}/{method}/{scenario}-{bkg_level}-{method}-{prefix}-result.fits.gz"
    output:
        "results/{scenario}/{bkg_level}/{prefix}/{method}/images/flux.png",
        "results/{scenario}/{bkg_level}/{prefix}/{method}/images/npred.png",
        "results/{scenario}/{bkg_level}/{prefix}/{method}/images/trace.png",
    log:
        "logs/plot-results/{scenario}-{bkg_level}-{method}-{prefix}.log"
    script:
        "../scripts/plot-results.py"