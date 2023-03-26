rule summarize_results:
    input:
        "config/{scenario}/{bkg_level}/{prefix}.yaml",
        "results/{scenario}/{bkg_level}/{prefix}/{method}/images/flux.png",
        "results/{scenario}/{bkg_level}/{prefix}/{method}/images/npred.png",
        "results/{scenario}/{bkg_level}/{prefix}/{method}/images/trace.png",
    output:
        "results/{scenario}/{bkg_level}/{prefix}/{method}/summary.rst"
    log:
        "logs/summarize-results/{scenario}-{bkg_level}-{method}-{prefix}.log"
    script:
        "../scripts/summarize-results.py"