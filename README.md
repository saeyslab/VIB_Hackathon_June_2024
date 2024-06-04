# VIB_Hackathon_June_2024
Code repository for the VIB Hackathon June 2024 on spatial omics. More info: https://hackmd.io/@berombau/BJetSxw8T

## Example working with SpatialData on Tier 1 HPC

Login to the hpc:

```bash
ssh {YOUR_USERNAME}@tier1.hpc.ugent.be
```

Set the following environment variables in `~/.bashrc`:

```
export SLURM_ACCOUNT='gpr_compute_starting_2024_011'
export SALLOC_ACCOUNT=$SLURM_ACCOUNT
export SBATCH_ACCOUNT=$SLURM_ACCOUNT
```

Start an interactive session, e.g.:

```bash
salloc --partition=debug_rome --nodes=1 --ntasks-per-node=4 --mem=16G --time=01:00:00
```

go to VSC_DATA_VO_USER if a VO is available, else go to $VSC_DATA_USER

```bash
cd $VSC_DATA_USER
```

Install [minconda](https://docs.anaconda.com/free/miniconda/#quick-command-line-install) there:

```bash
mkdir -p miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3/miniconda.sh
bash miniconda3/miniconda.sh -b -u -p miniconda3
rm -rf miniconda3/miniconda.sh
miniconda3/bin/conda init bash
```

Set `libmamba` as the solver:

```bash
conda install -n base conda-libmamba-solver
conda config --set solver libmamba
```

Create the environment:

```bash
conda env create -f env_python.yml
```

Activate the environment:

```bash
conda activate hackathon_multi_omics_2024
```

Try running the example [script](./examples/read_spatialdata_from_s3bucket.py):

```bash
python read_spatialdata_from_s3bucket.py
```

This script should create the file `test.png` in the folder `figures`.
