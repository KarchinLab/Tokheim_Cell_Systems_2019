# CHASMplus

With the ever-growing pace of DNA sequencing of human tumors, the total number of detected mutations in cancer continues to accelerate. However, only a few mutations in each tumor may actually “drive” the growth of cancer, some of which can have value for diagnostic, prognostic, or therapeutic purposes. Based on a new rigorous statistical analysis of The Cancer Genome Atlas (TCGA), we find a prominent emerging role for rare missense mutations predicted to be “drivers” of cancer, which may have potential implications for genome-driven precision oncology, since rare driver mutations that are putatively actionable could be newly observed in a patient, thus requiring personalized modeling and assessment. To extend beyond the TCGA, we provide a systematic resource to assess such newly observed missense mutations as cancer drivers. Lastly, we assess the driver landscape of human cancers and find that discovery for some cancer types are already approaching saturation.

In this Jupyter Notebook, we present a new statistically rigorous method, CHASMplus, for predicting the driver status of missense mutations. After careful benchmarking, we applied CHASMplus to 8,657 sequenced tumors from The Cancer Genome Atlas (TCGA) spanning 32 types of cancer. We explore the role for rare driver missense mutations in cancer and, when possible, relate predictions to supporting functional evidence. 

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

The notebooks use data and results available [here](http://karchinlab.org/data/CHASMplus/Tokheim_2018.tar.gz).
Place the data in the top-level directory of this repository.
The scores for each method in the benchmark are found in the `CHASMplus/data/benchmark` folder. If you are interested
in obtaining the full list of mutations used for training CHASMplus, it is available [here](http://karchinlab.org/data/CHASMplus/formatted_training_list.txt.gz).

## Citation

We will update the citation information upon publication.
