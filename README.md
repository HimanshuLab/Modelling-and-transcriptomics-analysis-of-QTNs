# Modelling-and-transcriptomics-analysis-of-QTNs
The entire repository is organized in two folders:
1) metabolic-modelling-of-SPO-QTNs
2) RNA_seq-analysis

### metabolic-modelling-of-SPO-QTNs
The folder **QTN_specific_metabolic_models** has all the QTN-specific metabolic models built in this study
__________________________________________________________________________

**Description of other files**

```QTN_specific_model_building.m```: This matlab code file builds all the QTN-specific models used in this study using the gene-expression data (**tpm_counts_Average.csv**) and the Yeast genome-scale model (**yeastGEM.mat**)

```optgp_sampling.ipynb```: This notebook is for sampling analyis and data visualisation



### RNA_seq-analysis
```Deseq2_analysis_up_down.R```: This is for differential gene expression analysis

```Rif_OOOO_TFs.R```: This is for regulatory impact factor analysis

### Acknowledgement
* [Centre for Integrative Biology and Systems medicinE](https://ibse.iitm.ac.in/)
* [Robert Bosch Centre for Data Science and Artificial Intelligence (RBCDSAI)](https://rbcdsai.iitm.ac.in/)
