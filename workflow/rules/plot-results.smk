
rule plot_results:
    input:
        unpack(get_datasets_filenames),
        filename_result="results/{scenario}/{bkg_level}/{prefix}/{method}/{scenario}-{bkg_level}-{method}-{prefix}-result.fits.gz",
        filename_flux_ref=get_filename_flux_ref,
        filename_config=get_filename_config,
    output:
        filename_image_flux="results/{scenario}/{bkg_level}/{prefix}/{method}/images/flux.png",
        filename_image_npred="results/{scenario}/{bkg_level}/{prefix}/{method}/images/npred.png",
        filename_image_trace="results/{scenario}/{bkg_level}/{prefix}/{method}/images/trace.png",
        filename_image_flux_thumbnail="results/{scenario}/{bkg_level}/{prefix}/{method}/images/flux-thumbnail.png",
        filename_metrics="results/{scenario}/{bkg_level}/{prefix}/{method}/metrics.yaml",
    log:
        "logs/plot-results/{scenario}-{bkg_level}-{method}-{prefix}.log"
    localrule: True
    script:
        "../scripts/plot-results.py"