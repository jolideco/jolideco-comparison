# Jolideco Comparison

Compare Jolideco against other deconvolution methods.

Get the data:

TODO: upload data to Zenodo and extend workflow to download it...or just submit...it is small.

Setup the environment:
```bash
mamba env create -f environment.yaml
```

Run comparison:
```bash
snakemake -c<n_process> --cluster workflow/cluster/qsub-submit.py --jobs <n_jobs>
```

Open webpage:
```
open site/index.html
```

Command to just rebuild the site:

```bash
snakemake -c8 -R site_render_index
```
