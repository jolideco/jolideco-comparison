rule site_render_index:
    input:
        expand("results/{scenario}/{bkg_level}/{prefix}/{method}/index.rst", scenario=config["scenarios"], bkg_level=config["bkg_levels"], prefix=config["prefixes"], method=config["methods"]),
        expand("results/{scenario}/{bkg_level}/{prefix}/images/counts.png", scenario=config["scenarios"], bkg_level=config["bkg_levels"], prefix=config["prefixes"]),
        expand("results/{scenario}/{bkg_level}/{prefix}/images/exposure.png", scenario=config["scenarios"], bkg_level=config["bkg_levels"], prefix=config["prefixes"]),
        expand("results/{scenario}/{bkg_level}/{prefix}/images/background.png", scenario=config["scenarios"], bkg_level=config["bkg_levels"], prefix=config["prefixes"]),
        expand("results/{scenario}/{bkg_level}/{prefix}/images/psf.png", scenario=config["scenarios"], bkg_level=config["bkg_levels"], prefix=config["prefixes"]),
    output:
        "results/index.rst"
    log:
        "logs/site-render-index.log"
    script:
        "../scripts/site-render-index.py"    