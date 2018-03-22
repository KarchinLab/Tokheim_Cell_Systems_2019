# CHASMplus

Sequencing studies have statistically implicated genetic drivers of human cancers by distinguishing these from the expected random accumulation of somatic mutations. However, prior work has coarsely focused on driver genes or regions, largely avoiding prediction of individual mutations. Here, we develop and rigorously validate CHASMplus to predict individual driver somatic missense mutations and show it exceeds state-of-the-art performance. Applied to 32 cancer types, CHASMplus identifies 3,527 unique drivers with four times higher prevalence of rare drivers than previously calculated. Our results indicate a complex relationship between the driver landscape of somatic missense mutations and each cancer type, some reveal a prominent role for rare drivers while others rely on common drivers and are already saturating discovery. The prevalence of rare cancer drivers has implications for future interpretation of cancer genomes and genome-driven oncology

## Jupyter notebooks

We have prepared our analysis into jupyter notebooks (.ipynb files). You can either view
them on github or execute the evaluation if you install [jupyter](http://jupyter.org/). We recommend that you first
view the [Introduction.ipynb](Introduction.ipynb) file for more details.

## Installation

We recommend you install all of the depdencies through conda. Please see the miniconda [installation page](https://conda.io/miniconda.html).

Next, install the dependencies needed to run these notebooks by the following commands:

```bash
$ conda env create -f environment.yml  # create environment for CHASMplus
$ source activate CHASMplus_jupyter  # activate environment for CHASMplus jupyter analysis
```

Remember to always activate the CHASMplus_jupyter environment in conda!

## Runing jupyter notebooks

Using the terminal, change directories to where you download this repository. Then start jupyter lab:

```bash
$ jupyter lab
```

## Data

The notebooks use data and results available [here](http://karchinlab.org/data/CHASMplus/Tokheim_NatGen_2018.tar.gz).
Place the data in the top-level directory of this repository.
The scores for each method in the benchmark are found in the `CHASMplus/data/benchmark` folder.

## Citation

We will update the citation information upon publication.
