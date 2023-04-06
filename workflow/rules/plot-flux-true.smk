rule plot_flux_true:
    input:
        get_filename_flux_ref,
        filename_config="config/{scenario}/bg1/chandra.yaml",
    output:
        filename_image_flux_ref="results/{scenario}/images/flux-true-thumbnail.png",
    log:
        "logs/plot-flux-true/{scenario}.log"
    localrule: True
    script:
        "../scripts/plot-flux-true.py"