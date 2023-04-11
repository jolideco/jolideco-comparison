rule plot_datasets:
    input:
        unpack(get_datasets_filenames),
        config="config/{scenario}/{bkg_level}/{prefix}.yaml",
    output:
        filename_image_counts="results/{scenario}/{bkg_level}/{prefix}/images/counts.png",
        filename_image_exposure="results/{scenario}/{bkg_level}/{prefix}/images/exposure.png",
        filename_image_background="results/{scenario}/{bkg_level}/{prefix}/images/background.png",
        filename_image_psf="results/{scenario}/{bkg_level}/{prefix}/images/psf.png",
        filename_image_counts_thumbnail="results/{scenario}/{bkg_level}/{prefix}/images/counts-thumbnail.png",
    log:
        "logs/plot-datasets/{scenario}-{bkg_level}-{prefix}.log"
    localrule: True
    script:
        "../scripts/plot-datasets.py"