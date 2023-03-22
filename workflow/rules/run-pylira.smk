rule run_jolideco:
    input:
        "table.txt"
    output:
        "plots/myplot.pdf"
    conda:
        "envs/jolideco.yaml"
    script:
        "scripts/run-jolideco.py"