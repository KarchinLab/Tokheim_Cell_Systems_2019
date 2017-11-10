# CHASM2

Sequencing studies have statistically implicated genetic drivers of human cancers by distinguishing these from the expected random accumulation of somatic mutations. However, prior work has coarsely focused on driver genes or regions, largely avoiding prediction of individual mutations. Here, we develop and rigorously validate CHASM2 to predict individual driver somatic missense mutations and show it exceeds state-of-the-art performance. Applied to 32 cancer types, CHASM2 identifies 3,527 unique drivers with four times higher prevalence of rare drivers than previously calculated. Our results indicate a complex relationship between the driver landscape of somatic missense mutations and each cancer type, some reveal a prominent role for rare drivers while others rely on common drivers and are already saturating discovery. We show experimentally that CHASM2 discriminates radiosensitivity mutations within the ATM gene, a phenotype of potential clinical relevance. The prevalence of rare cancer drivers has implications for future interpretation of cancer genomes and genome-driven oncology

## Jupyter notebooks

We have prepared our analysis into jupyter notebooks (.ipynb files). You can either view
them on github or execute the evaluation if you install [jupyter](http://jupyter.org/). We recommend that you first
view the [Introduction.ipynb](LINK to come) file for more details.

## Installation

We recommend you install all of the depdencies through conda. Please see the miniconda [installation page](https://conda.io/miniconda.html).

Next, install the dependencies needed to run these notebooks by the following commands:

```bash
$ conda env create -f environment.yml  # create environment for CHASM2
$ source activate CHASM2_jupyter  # activate environment for CHASM2 jupyter analysis
```

Remember to always activate the CHASM2_jupyter environment in conda!

## Runing jupyter notebooks

Using the terminal, change directories to where you download this repository. Then start jupyter lab:

```bash
$ jupyter lab
```

## Citation

You will update the citation information upon publication.
