# PertInInt Processing Steps

Welcome to the internal processing steps for PertInInt! Here, you will find instructions to recreate all input files required by PertInInt from scratch. If you find code from this repository useful, please cite:

> Kobren, S.N., Chazelle, B. and Singh, M. (2019). "An integrative approach to identify preferentially altered interactions in human cancers." *Manuscript in preparation.*

## Configuration

**Before starting,** you *must* first create a configuration file called **config.py** in this git repository. You might run this code on your local work machine, laptop, cluster, etc., and where you store your data files will be different in each location. For instance, on the Princeton cluster, my `config.py` file looks like:

```python
data_path = '/Genomics/grid/users/snadimpa/data/'
bin_dir = '/Genomics/grid/users/snadimpa/bin/'
```

## Critical Steps in the PertInInt pipeline

**The following Wiki pages describe all major steps in the PertInInt pipeline.** If you have any questions, please contact Shilpa N. Kobren at snadimpa@alumni.princeton.edu. 

1. [Verify protein sequences from Ensembl](https://github.com/Singh-Lab/pertinint-internal/wiki/Verify-protein-sequences-from-Ensembl)
2. [Define functional tracks](https://github.com/Singh-Lab/pertinint-internal/wiki/Define-functional-tracks)
   - [interaction tracks](https://github.com/Singh-Lab/pertinint-internal/wiki/Interaction-tracks)
   - [domain tracks](https://github.com/Singh-Lab/pertinint-internal/wiki/Domain-tracks)
   - [conservation tracks](https://github.com/Singh-Lab/pertinint-internal/wiki/Conservation-tracks)
   - [natural variation tracks](https://github.com/Singh-Lab/pertinint-internal/wiki/Natural-variation-tracks)
3. [Download TCGA cancer data](https://github.com/Singh-Lab/pertinint-internal/wiki/Download-TCGA-cancer-data)
4. [Preprocess functional tracks for use by PertInInt](https://github.com/Singh-Lab/pertinint-internal/wiki/Preprocess-functional-tracks)
5. [Create files to annotate output](https://github.com/Singh-Lab/pertinint-internal/wiki/Annotate-gene-names-and-drivers)
6. [Run PertInInt](https://github.com/Singh-Lab/PertInInt)
