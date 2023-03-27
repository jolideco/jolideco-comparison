rule prepare_site:
    input:
       expand("results/{scenario}/{bkg_level}/{prefix}/{method}/index.rst", scenario=config["scenarios"], bkg_level=config["bkg_levels"], prefix=config["prefixes"], method=config["methods"]),
    output:
        "results/index.rst"
    log:
        "logs/prepare-site.log"
    script:
        "../scripts/prepare-site.py"    