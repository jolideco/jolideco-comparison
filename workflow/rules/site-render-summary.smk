rule site_render_summary:
    input:
        filename_config="config/{scenario}/{bkg_level}/{prefix}.yaml",
        filename_image_flux="results/{scenario}/{bkg_level}/{prefix}/{method}/images/flux.png",
        filename_image_npred="results/{scenario}/{bkg_level}/{prefix}/{method}/images/npred.png",
        filename_image_trace="results/{scenario}/{bkg_level}/{prefix}/{method}/images/trace.png",
        filename_result="results/{scenario}/{bkg_level}/{prefix}/{method}/{scenario}-{bkg_level}-{method}-{prefix}-result.fits.gz",
        filename_metrics="results/{scenario}/{bkg_level}/{prefix}/{method}/metrics.yaml",
    output:
        filename_index="results/{scenario}/{bkg_level}/{prefix}/{method}/index.rst",
    log:
        "logs/site-render-summary/{scenario}-{bkg_level}-{method}-{prefix}.log"
    localrule: True
    script:
        "../scripts/site-render-summary.py"