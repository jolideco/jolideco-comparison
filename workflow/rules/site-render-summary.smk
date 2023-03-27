rule site_render_summary:
    input:
        "config/{scenario}/{bkg_level}/{prefix}.yaml",
        "results/{scenario}/{bkg_level}/{prefix}/{method}/images/flux.png",
        "results/{scenario}/{bkg_level}/{prefix}/{method}/images/npred.png",
        "results/{scenario}/{bkg_level}/{prefix}/{method}/images/trace.png",
    output:
        "results/{scenario}/{bkg_level}/{prefix}/{method}/index.rst",
    log:
        "logs/site-render-summary/{scenario}-{bkg_level}-{method}-{prefix}.log"
    script:
        "../scripts/site-render-summary.py"