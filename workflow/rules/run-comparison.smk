rule run_comparison:
    input:
        "config/{name}/{bkg_level}/{prefix}.yaml"
    output:
        "results/{name}/{bkg_level}/{prefix}/{method}/{name}-{bkg_level}-{method}-{prefix}-result.fits.gz"
    log:
        "logs/run-comparison/{name}-{bkg_level}-{method}-{prefix}.log"
    script:
        "../scripts/run-comparison.py"