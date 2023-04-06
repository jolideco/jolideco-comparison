rule run_comparison:
    input:
        unpack(get_datasets_filenames),
        config="config/{scenario}/{bkg_level}/{prefix}.yaml",
    output:
        "results/{scenario}/{bkg_level}/{prefix}/{method}/{scenario}-{bkg_level}-{method}-{prefix}-result.fits.gz"
    log:
        "logs/run-comparison/{scenario}-{bkg_level}-{method}-{prefix}.log"
    script:
        "../scripts/run-comparison.py"