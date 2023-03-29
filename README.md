# bacterial shotgun assembly  

======================================================

Snakemake workflow for de novo assembly of bacterial shotgun seq. Please not that this workflow does not perform additional investigation of contigs and scaffolds as well as manual curation of the assembled genome. 

======================================================

run on your local machine:  
1) set up and configure units.tsv in config/ folder  
2) add your genome & annotation of choice to the resources/ folder  
3) identify the run folder (where the sequencer saved your run data)  
4) do a dry-run and check if the right jobs are to be run  

```
snakemake -npr --config runID="20230125_RNASeq"

```

5) run the job locally (choose the correct runID - this can be found in the SampleSheet.csv)
```
time snakemake -j48 --config runID="20230125_RNASeq"
```

======================================================

## conda/mamba and other [dependencies](https://github.com/schimar/rna_fusion_quant/blob/main/workflow/envs/s7.yaml)   

create environment from yaml file (in envs/):
```
# create the environment (note that conda & mamba have to be installed for this to work):

mamba env create -f envs/s6.yaml

# with this, you can activate the environment with all [dependencies](https://github.com/schimar/rna_fusion_quant/blob/main/workflow/envs/s7.yaml):
conda activate smk6


# if you've added new software to install to the conda environment, then you can update:
mamba env update -f envs/s7.yaml
```

## with SLURM, submit the job:
```
sbatch code/clusterSnakemake.sh
```


### Deprecated (or rather: not used here): submit on PBS
```
qsub code/clusterSnakemake.pbs
```




