rule site_build:
    input:
        "results/index.rst"
    output:
        directory("site/")
    log:
        "logs/site-build/sphinx.log"
    localrule: True
    shell:
        "cp workflow/site/conf.py results/;"
        "cp -r workflow/site/_static results/;"
        "sphinx-build results site;"