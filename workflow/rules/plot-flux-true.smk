rule plot_flux_true:
    input:
        "config/{scenario}/bg1/chandra.yaml"
    output:
        "results/{scenario}/images/flux-true-thumbnail.png",
    log:
        "logs/plot-flux-true/{scenario}.log"
    localrule: True
    script:
        "../scripts/plot-flux-true.py"