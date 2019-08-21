# Snakemake workflow: 201904_MPESeq

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.3.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/201904_MPESeq.svg?branch=master)](https://travis-ci.org/snakemake-workflows/201904_MPESeq)


## Authors

* Benjamin Fair (@bfairkun)

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, clone it:
```
git clone https://github.com/bfairkun/201907_mpe-seq_processing.git
```
If you intend to modify and further develop this workflow, fork this repository. Please consider providing any generally applicable modifications via a pull request.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml` and the `samples.tsv` files. Configure cluster settings in `cluster-config.json`. If you want to run snakemake from a compute node using the provided `snakemake.sbatch` script, edit it accordingly (specifically the account for submitting slurm jobs unless).

### Step 3: Prepare the environments

Most of necessary software, including snakemake, can be installed with conda using the `envs/general_environment.yaml` file. There may exist additional environments that snakemake may automatically switch to on a rule-by-rule basis that need to be installed.

 ```
 # create the general environment called rna_seq_processing from the provided yaml file
 conda env create -f envs/general_environment.yaml
 
 # activate the environment
 conda activate rna_seq_processing
 
 # create the conda environments specified in the workflow rules using snakemake
 snakemake --use-conda --create-envs-only
 ```

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --cluster --cluster-config cluster-config.json --cluster "sbatch --partition={cluster.partition} --job-name={cluster.name} --output=/dev/null --job-name={cluster.name} --nodes={cluster.n} --mem={cluster.mem}"

or by executing the included sbatch script to execute the snakemake process from a cluster

    sbatch snakemake.sbatch

See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.
