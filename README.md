# VIB_Hackathon_June_2024
Code repository for the VIB Hackathon June 2024 on spatial omics. More info: https://hackmd.io/@berombau/BJetSxw8T.

## Compile pdf

See instructions at https://github.com/biohackrxiv/bhxiv-gen-pdf.

Update the `paper.pdf` using the `paper.md`.
```
docker run --rm -it -v $(pwd):/work -w /work biohackrxiv/gen-pdf:local gen-pdf .
```

## Example working with SpatialData on Tier 1 HPC

Login to the hpc:

```bash
ssh {YOUR_USERNAME}@tier1.hpc.ugent.be
```

Set the following environment variables in `~/.bashrc`:

```
export SLURM_ACCOUNT='starting_2024_011'
export SALLOC_ACCOUNT=$SLURM_ACCOUNT
export SBATCH_ACCOUNT=$SLURM_ACCOUNT
```

Start an interactive session, e.g.:

```bash
qsub -I -l nodes=1:ppn=16,mem=4G
```

go to `VSC_DATA_VO_USER` if a VO is available, else go to `VSC_DATA_USER`

```bash
cd $VSC_DATA_USER
```

Clone the repository

```bash
git clone https://github.com/saeyslab/VIB_Hackathon_June_2024.git
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
cd $VSC_DATA_USER/VIB_Hackathon_June_2024
conda env create -f env_python.yml
```

Activate the environment:

```bash
conda activate hackathon_multi_omics_2024
```

Try running the example [script](./examples/read_spatialdata_from_s3bucket.py):

```bash
python examples/read_spatialdata_from_s3bucket.py
```

This script should create the file `test.png` in the current working directory.

We also include a simple `.pbs` [script](./examples/run_hackathon.pbs), that can be submitted to the job queue, e.g.:

```bash
module swap cluster/dodrio/cpu_rome
qsub examples/run_hackathon.pbs
```

This will also create the file `test.png` in the current working directory.

## Example working with hydra on Tier 1 HPC

We provide an example on how to launch jobs on the HPC using [hydra](https://hydra.cc/docs/intro/) and [hydra-submitit](https://hydra.cc/docs/plugins/submitit_launcher/). After cloning the repository and installing the conda environment (see previous section) do the following to launch a job:

```bash
module swap cluster/dodrio/cpu_rome
cd $VSC_DATA_USER/VIB_Hackathon_June_2024/examples
chmod +x run_hydra.sh
./run_hydra.sh
```

A job will be queued, and in the folder `examples` a directory `multirun` will be created, containing the logs of your job.

You can play with this example locally if you change the following section of [`default.yaml`](./examples/configs/default.yaml) from this

```yaml
defaults:
  - override /hydra/launcher: submitit_slurm
```

to this

```yaml
defaults:
  - override /hydra/launcher: submitit_local
```

However, local testing of `hydra-submitit` will not work on Windows.

For more complex examples on how to work with `hydra`, we refer to [here](https://github.com/saeyslab/hydra_hpc_example) and [here](https://github.com/saeyslab/napari-sparrow).

## Steps to Connect to Tier 1 HPC with VS Code Remote - SSH


### Install VS Code and the Remote - SSH Extension:

    - Open VS Code.

    - Go to the Extensions view by clicking on the Extensions icon in the Activity Bar on the side of the window or pressing Ctrl+Shift+X.

    - Search for "Remote - SSH" and install it.

### Configure SSH in VS Code:

    - Press F1 to open the Command Palette in VS Code.

    - Type Remote-SSH: Open Configuration File... and select it.

    - Choose the SSH configuration file to edit. It’s usually located at ~/.ssh/config.

    - Add the following configuration to the file:

```yaml
Host hpc_tier1
    HostName tier1.hpc.ugent.be
    User {YOUR_USERNAME}
    ForwardAgent yes
```

### Connect to the Remote Server:

    - Press F1 to open the Command Palette.

    - Type Remote-SSH: Connect to Host... and select it.

    - You should see hpc_tier1 in the list. Select it.

    - VS Code will open a new window and start connecting to the remote server.
    
    - You might be prompted for your SSH key passphrase or password.