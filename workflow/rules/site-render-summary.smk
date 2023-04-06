rule site_render_summary:
    input:
        filename_config="config/{scenario}/{bkg_level}/{prefix}.yaml",
        filename_image_flux="results/{scenario}/{bkg_level}/{prefix}/{method}/images/flux.png",
        "results/{scenario}/{bkg_level}/{prefix}/{method}/images/npred.png",
        "results/{scenario}/{bkg_level}/{prefix}/{method}/images/trace.png",
        "results/{scenario}/{bkg_level}/{prefix}/{method}/{scenario}-{bkg_level}-{method}-{prefix}-result.fits.gz",
        "results/{scenario}/{bkg_level}/{prefix}/{method}/metrics.yaml",
    output:
        "results/{scenario}/{bkg_level}/{prefix}/{method}/index.rst",
    log:
        "logs/site-render-summary/{scenario}-{bkg_level}-{method}-{prefix}.log"
    localrule: True
    script:
        "../scripts/site-render-summary.py"