def get_flux_true_filename(wildcards):
    """"Get flux ground truth filename"""
    # grpound truth is the same for Chandra and XMM and all bkg levels
    filename = get_filenames("chandra", bkg_level="", name=wildcards.scenario, quantity="flux")[0]
    return filename

rule plot_flux_true:
    input:
        get_flux_true_filename,
        config="config/{scenario}/bg1/chandra.yaml",
    output:
        "results/{scenario}/images/flux-true-thumbnail.png",
    log:
        "logs/plot-flux-true/{scenario}.log"
    localrule: True
    script:
        "../scripts/plot-flux-true.py"