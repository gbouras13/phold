# Installation

The only way to install `phold` is from source for now. 

Pypi and (hopefully) conda installations will be available soon. 

The only required non-Python dependency is [Foldseek](https://github.com/steineggerlab/foldseek). To install `phold` in a conda environment using [mamba](https://github.com/conda-forge/miniforge):

```
mamba create -n pholdENV pip foldseek python=3.11
conda activate pholdENV
git clone https://github.com/gbouras13/phold.git
cd phold 
pip install -e .
```

## Torch 

To utilise `phold` with GPU, a GPU compatible version of `pytorch` must be installed. 

If it is not automatically installed via the pip/conda installation, please see [this link](https://pytorch.org) for more instructions on how to install `pytorch`. 

If you have an older version of the CUDA driver installed on your NVIDIA GPU, then you might find [this link useful](https://pytorch.org/get-started/previous-versions/).

Phold has been tested on NVIDIA GPUs (A100, RTX4090) and AMD GPUs (Radeon). 

Installation on AMD GPUs requires `torch` compatible with rocm e.g.

```
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/rocm5.7
```

# Database Installation

To download and install the `phold` database

```
phold install
```

If you would like to specify a particular location for the database, please use `-d`

```
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

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

We would recommend installing `phold` into a fresh environment. Assuming you installed miniforge, to create a environment called `pholdENV` with `phold` installed:

* To create a conda environment called `pholdENV` with foldseek installed

```
conda create -n pholdENV foldseek pip
```

* To activate the environment

```
conda activate pholdENV
```

* To install `phold`

```
pip install phold
```

* Once that has finished downloading and installing, you can check installation worked using:

```
phold -h
phold install
```