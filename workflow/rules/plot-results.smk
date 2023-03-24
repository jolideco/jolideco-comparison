rule plot_results:
    input:
        "results/{name}/{bkg_level}/{prefix}/{method}/{name}-{bkg_level}-{method}-{prefix}-result.fits.gz"
    output:
        "results/{name}/{bkg_level}/{prefix}/{method}/images/flux.png",
        "results/{name}/{bkg_level}/{prefix}/{method}/images/npred.png",
        "results/{name}/{bkg_level}/{prefix}/{method}/images/trace.png",
    log:
        "logs/plot-results/{name}-{bkg_level}-{method}-{prefix}.log"
    script:
        "../scripts/plot-results.py"