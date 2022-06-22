
To run this workflow:

```
bsub nextflow -C nf_config/nextflow.config \
              run workflows/run_protein_prediction.nf \
              -profile lsf -resume
```
