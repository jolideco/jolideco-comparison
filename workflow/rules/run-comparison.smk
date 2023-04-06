rule run_comparison:
    input:
        unpack(get_datasets_filenames),
        config="config/{scenario}/{bkg_level}/{prefix}.yaml",
    output:
        expand("results/{{scenario}}/{{bkg_level}}/{{prefix}}/{method}/{{scenario}}-{{bkg_level}}-{method}-{{prefix}}-result.fits.gz", method=config["methods"])
    log:
        expand("logs/run-comparison/{{scenario}}-{{bkg_level}}-{method}-{{prefix}}.log", method=config["methods"])
    script:
        "../scripts/run-comparison.py"