#!/usr/bin/env python3
import os
import sys
from pathlib import Path

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

os.system(f"cp {jobscript} print.sh")

print(jobscript)
print(job_properties["jsdhgdh"])

# os.system("qsub {script}".format(script=jobscript))
