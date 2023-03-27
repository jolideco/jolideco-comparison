rule build_site:
    input:
        "results/index.rst"
    output:
        directory("site/")
    log:
        "logs/build-site/sphinx.log"
    shell:
        "sphinx-build results site"       