rule site_build:
    input:
        "results/index.rst"
    output:
        "site/index.html"
    log:
        "logs/site-build/sphinx.log"
    localrule: True
    shell:
        "cp workflow/site/conf.py results/;"
        "cp -r workflow/site/_static results/;"
        "sphinx-build results site;"