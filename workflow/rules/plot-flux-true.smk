rule plot_flux_true:
    output:
        "results/{scenario}/images/flux-true-thumbnail.png",
    log:
        "logs/plot-flux-true/{scenario}.log"
    script:
        "../scripts/plot-flux-true.py"