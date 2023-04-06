rule plot_datasets:
    input:
        unpack(get_datasets_filenames),
        config="config/{scenario}/{bkg_level}/{prefix}.yaml",
    output:
        "results/{scenario}/{bkg_level}/{prefix}/images/counts.png",
        "results/{scenario}/{bkg_level}/{prefix}/images/exposure.png",
        "results/{scenario}/{bkg_level}/{prefix}/images/background.png",
        "results/{scenario}/{bkg_level}/{prefix}/images/psf.png",
    log:
        "logs/plot-datasets/{scenario}-{bkg_level}-{prefix}.log"
    localrule: True
    script:
        "../scripts/plot-datasets.py"