# Installation

The best way to install `phold` is using [mamba](https://github.com/conda-forge/miniforge), as this will install [Foldseek](https://github.com/steineggerlab/foldseek) (the only non-Python dependency) along with the Python dependencies.

To install `phold` using [mamba](https://github.com/conda-forge/miniforge):

```bash
mamba create -n pholdENV -c conda-forge -c bioconda phold 
```

To utilise `phold` with GPU, a GPU compatible version of `pytorch` must be installed. By default conda/mamba will install a CPU-only version. 

Therefore, if you have an NVIDIA GPU, please try:

```bash
mamba create -n pholdENV -c conda-forge -c bioconda phold pytorch=*=cuda*
```

## Pip

You can also install the `phold` using pip.

```bash
pip install phold
```

You will need to have [Foldseek](https://github.com/steineggerlab/foldseek) installed and available in the $PATH.

## Source

You can install the latest version of `phold` with potentially untested and unreleased changes into a conda environment using [mamba](https://github.com/conda-forge/miniforge) as follows:

```bash
mamba create -n pholdENV pip foldseek python=3.11
conda activate pholdENV
git clone https://github.com/gbouras13/phold.git
cd phold 
pip install -e .
```

## Torch 

To utilise `phold` with GPU, a GPU compatible version of `pytorch` must be installed. 

If it is not automatically installed via the installation methods above, please see [this link](https://pytorch.org) for more instructions on how to install `pytorch`. 

If you have an older version of the CUDA driver installed on your NVIDIA GPU, then you might find [this link useful](https://pytorch.org/get-started/previous-versions/).

Phold has been tested on NVIDIA GPUs (A100, RTX4090) and AMD GPUs (Radeon). 

Installation on AMD GPUs requires a version of `torch` compatible with rocm e.g.

```bash
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/rocm5.7
```

# Database Installation

To download and install the `phold` database

```bash
phold install
```

If you would like to specify a particular location for the database (e.g. if you use `phold` on a shared server), please use `-d`

```bash
phold install -d <path/to/databse_dir>
```

* Note: You will need at least 8GB of free space (the `phold` databases including ProstT5 are 7.7GB uncompressed).

# Beginner Conda Installation

If you are new to using the command-line, please install conda using the following instructions.

First install some flavour of [Anaconda](https://www.anaconda.com/products/distribution). 

There are lots of options but the best in our opinion is miniforge as this will automatically install mamba, which is much faster than base conda:

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
mamba create -n pholdENV -c conda-forge -c bioconda phold pytorch=*=cuda*
```

If you have a Mac that runs Apple Silicon (M1/M2/M3),  please try:

```bash
mamba create -n pholdENV python==3.11  
conda activate pholdENV
mamba install pytorch::pytorch torchvision torchaudio -c pytorch 
mamba install -c conda-forge -c bioconda phold 
```

If you don't have a GPU:

```bash
mamba create -n pholdENV -c conda-forge -c bioconda phold 
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