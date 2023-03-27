rule site_build:
    input:
        "results/index.rst"
    output:
        directory("site/")
    log:
        "logs/site-build/sphinx.log"
    shell:
        "sphinx-build results site"       