# Installation

## Conda

The easiest way to install `phold` is via conda/mamba. This will install all python dependencies and [Foldseek](https://github.com/steineggerlab/foldseek), the only non-python dependency of `phold` all at once.

```
maba install -c bioconda phold
```

## Pip

You can also install `phold` using pip

```
pip install phold
```

You will still need to install [Foldseek](https://github.com/steineggerlab/foldseek) separately.

## Source

Alternatively, the development version of `phold` (which may include new, untested features) can be installed manually via github. 

```
git clone https://github.com/gbouras13/pharokka.git
cd phold
pip install -e .
phold --help
```

You will still need to install [Foldseek](https://github.com/steineggerlab/foldseek) separately.


# Database Installation

To download and install the `phold` database

```
phold install -d <path to database directory>
```
<!-- 
If you would like to specify a different database directory (recommended), that can be achieved as follows:

`install_databases.py -o <path/to/databse_dir>`

If this does not work, you an alternatively download the databases from Zenodo at https://zenodo.org/record/8276347/files/pharokka_v1.4.0_databases.tar.gz and untar the directory in a location of your choice.

If you prefer to use the command line:

```
wget "https://zenodo.org/record/8267900/files/pharokka_v1.4.0_databases.tar.gz"
tar -xzf pharokka_v1.4.0_databases.tar.gz
```

which will create a directory called "pharokka_v1.4.0_databases" containing the databases. -->

Beginner Conda Installation
----

If you are new to using the command-line, please install conda using the following instructions.

First install [Anaconda](https://www.anaconda.com/products/distribution). There are lots of options but the best in our opinion is miniforge as this will automatically install mamba, which is much faster than base conda:

   * [miniforge](https://github.com/conda-forge/miniforge).
  
Please follow the instructions at the links to install based on your computer architecture. 

After your installation is complete, you should add the following channels to your conda configuration:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

We would recommend installing `phold` into a fresh environment. Assuming you installed miniforge, to create a environment called `pholdENV` with `phold` installed:

* To create a conda environment called `pholdENV`

```
conda create -n pholdENV
```

* To activate the environment

```
conda activate pholdENV
```

* To install phold

```
mamba install -c bioconda phold
```

* Once that has finished downloading and installing, you can check installation worked using:

```
phold -h
```