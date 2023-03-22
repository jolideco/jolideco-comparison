rule run_comparison:
    input:
        "config/{name}/{bkg_level}/{prefix}.yaml"
    output:
        directory("results/{name}/{bkg_level}/{method}")
    log:
        "logs/{name}/{bkg_level}/{method}.log"
    script:
        "scripts/run-comparison.py"