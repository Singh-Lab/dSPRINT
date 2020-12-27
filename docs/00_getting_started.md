## Getting Started

dSPRINT can be run on a any x64 Linux machine (either a local desktop/laptop or on a cluster). The only requirements
are a 64 bit [Anaconda installation](https://www.anaconda.com/products/individual), and plenty of disk space
(around 400G in a typical case).

### Install Anaconda

If you don't already have Anaconda installed, use the instructions at
[https://www.anaconda.com/products/individual](https://www.anaconda.com/products/individual) to download and install
Anaconda for your platform. (Note: The lite version of Anaconda, [miniconda](https://docs.conda.io/en/latest/miniconda.html) will work just fine).

On a cluster environment, you may already have Anaconda available by way of 
[Environment Modules](https://modules.readthedocs.io/en/latest/). If you run:

```
module avail anaconda
```

or 

```
module avail conda
```

and you see one or more entries returned, you can execute `module load anaconda` or `module load conda` to bring the
`conda` command in your PATH. If this is the first time you're using `conda`, you will also want to execute `conda init`,
which is a required step for recent versions of conda.


### Set up environment

After cloning the code in this repository:

```
git clone https://github.com/vineetbansal/dsprint-pipeline.git
```

`cd` to the folder where you clone the repository (this folder will have the file `environment.yml` and `config.json`,
 among others), and create a new conda environment where you can run the dSPRINT pipeline:

```
cd dsprint-pipeline
conda env create -f environment.yml
```

Activate the newly created environment, which is called `dsprint3`:

```
conda activate dsprint
```

Proceed to the [Downloading Data](01_download_data.md) tutorial.
