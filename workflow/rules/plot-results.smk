rule plot_results:
    input:
        "results/{scenario}/{bkg_level}/{prefix}/{method}/{scenario}-{bkg_level}-{method}-{prefix}-result.fits.gz"
    output:
        "results/{scenario}/{bkg_level}/{prefix}/{method}/images/flux.png",
        "results/{scenario}/{bkg_level}/{prefix}/{method}/images/npred.png",
        "results/{scenario}/{bkg_level}/{prefix}/{method}/images/trace.png",
        "results/{scenario}/{bkg_level}/{prefix}/{method}/images/flux-thumbnail.png",
        "results/{scenario}/{bkg_level}/{prefix}/{method}/metrics.yaml",
    log:
        "logs/plot-results/{scenario}-{bkg_level}-{method}-{prefix}.log"
    localrule: True
    script:
        "../scripts/plot-results.py"