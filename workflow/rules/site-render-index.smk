rule site_render_index:
    input:
       expand("results/{scenario}/{bkg_level}/{prefix}/{method}/index.rst", scenario=config["scenarios"], bkg_level=config["bkg_levels"], prefix=config["prefixes"], method=config["methods"]),
    output:
        "results/index.rst"
    log:
        "logs/site-render-index.log"
    script:
        "../scripts/site-render-index.py"    