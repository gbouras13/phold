# Installation

The best way to install `phold` is using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html), as this will install [Foldseek](https://github.com/steineggerlab/foldseek) (the only non-Python dependency) along with the Python dependencies.

We would highly recommend installing conda via [miniforge](https://github.com/conda-forge/miniforge).

To install `phold`:

```bash
conda create -n pholdENV -c conda-forge -c bioconda phold 
```

To utilise `phold` with GPU, a GPU compatible version of `pytorch` must be installed. By default conda will install a CPU-only version. 

Therefore, if you have an NVIDIA GPU, please try:

```bash
conda create -n pholdENV -c conda-forge -c bioconda phold pytorch=*=cuda*
```

## Pip

You can also install `phold` using pip.

```bash
pip install phold
```

You will need to have [Foldseek](https://github.com/steineggerlab/foldseek) (ideally v10.941cd33) installed and available in the $PATH.

## Source

You can install the latest version of `phold` with potentially untested and unreleased changes into a conda environment using [conda](https://github.com/conda-forge/miniforge) as follows:

```bash
conda create -n pholdENV pip foldseek python=3.13
conda activate pholdENV
git clone https://github.com/gbouras13/phold.git
cd phold 
pip install -e .
```

## Mac (M1/M2/M3)

If you have a Mac that runs Apple Silicon (M1/M2/M3), please try:

```bash
mamba create -n pholdENV python==3.10
conda activate pholdENV
mamba install pytorch::pytorch torchvision torchaudio -c pytorch 
mamba install -c conda-forge -c bioconda phold 
```

## pytorch 

To utilise `phold` with GPU, a GPU compatible version of `pytorch` must be installed. 

If it is not automatically installed via the installation methods above, please see [this link](https://pytorch.org) for more instructions on how to install `pytorch`. 

If you have an older version of the CUDA driver installed on your NVIDIA GPU, then you might find [this link useful](https://pytorch.org/get-started/previous-versions/).

Phold has been tested on NVIDIA (A100, RTX4090), AMD (MI250) and Mac (M1 Pro) GPUs. 

Installation on AMD GPUs requires a version of `torch` compatible with rocm e.g.

```bash
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/rocm5.7
```

# Database Installation

To download and install the `phold` database. Database downloads are multi-threaded (as of v1.0.0).

```bash
phold install -t <threads>
```

If you would like to specify a particular location for the database (e.g. if you use `phold` on a shared server), please use `-d`

```bash
phold install -d <path/to/databse_dir> -t <threads>
```

* Note: You will need at least 8GB of free space (the `phold` databases including ProstT5 are 7.7GB uncompressed).

If you have an NVIDIA GPU available, you may wish to accelerate Foldseek using GPU. To do this, you will need to format the databases appropriately as follows

```bash
phold install -d <path/to/databse_dir> --foldseek_gpu -t <threads>
```

If you would like to install the Phold DB 3.16M (i.e. including 1.8M extra enVhog and efam proteins without PHROG assignment or functional annotation label), you can use `--extended_db`

```bash
Usage: phold install [OPTIONS]

  Installs ProstT5 model and phold database

Options:
  -h, --help             Show this message and exit.
  -V, --version          Show the version and exit.
  -d, --database TEXT    Specific path to install the phold database
  --foldseek_gpu         Use this to enable compatibility with Foldseek-GPU
                         acceleration
  --extended_db          Download the extended Phold DB 3.16M including 1.8M
                         efam and enVhog proteins without functional labels
                         instead of the default Phold Search 1.36M. Using the
                         extended database will likely marginally reduce
                         functional annotation sensitivity and increase
                         runtime, but may find more hits overall i.e.
                         including to efam and enVhog proteins that have no
                         functional labels.
  -t, --threads INTEGER  Number of threads  [default: 1]
```

# Beginner Conda Installation

If you are new to using the command-line, please install conda using the following instructions.

First install some flavour of [Anaconda](https://www.anaconda.com/products/distribution). 

There are lots of options but the best in our opinion is miniforge as this will automatically use the mamba solver, which has now been integrated into conda:

   * [miniforge](https://github.com/conda-forge/miniforge).
  
Please follow the instructions at the links to install based on your computer architecture. 

After your installation is complete, you should add the following channels to your conda configuration:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

We would recommend installing `phold` into a fresh environment. Assuming you installed miniforge, to create an environment called `pholdENV` with `phold` installed (assuming you have an NVIDIA GPU):

```bash
conda create -n pholdENV -c conda-forge -c bioconda phold pytorch=*=cuda*
```
If you don't have a GPU:

```bash
conda create -n pholdENV -c conda-forge -c bioconda phold 
```

Then activate the environment

```bash
conda activate pholdENV
```

You can then check installation worked and download the `phold` databases:

```bash
phold -h
phold install
```

See the [tutorial](https://phold.readthedocs.io/en/latest/tutorial/) for more information on how to run `phold`.