#!/usr/bin/env python3
import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

name_template = "{scenario}-{bkg_level}-{prefix}"
name = job_properties["rule"] + name_template.format(**job_properties["wildcards"])

ARGS = {
    "-S": "/bin/bash",
    "-pe": "mthread 8",
    "-q": "sThC.q",
    "-l": "mres=64G,h_data=8G,h_vmem=8G",
    "-j": "y",
    "-N": f"{name}",
    "-M": "axel.donath@cfa.harvard.edu",
}

args = " ".join(f"{key} {value}" for key, value in ARGS.items())

os.system(f"qsub {args} {jobscript}")
