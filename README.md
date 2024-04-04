# Interpretable Machine Learning Decodes Soil Microbiome’s Response to Drought Stress

This repository contains the code used for 16S rRNA-based metagenomics analysis and ML applications for drought stress identification in the soil metagenome using a dataset from Naylor *et al.* (2017) further refered to as the **grass-drought** dataset.

- `SoilMicrobiomeDroughtML/`
  - `Metagenomics_Analysis.md`: Markdown showing how Data Processing, Diversity Analysis and Differential Abundance Analysis was performed on the grass-drought dataset
  - `Metagenomics_Analysis_files/`: Folder containing plots to be displayed in `Metagenomics_Analysis.md`
  - `Marker_Taxa.md`: Markdown showing the comparison of significant taxa from DAA and important taxa from interpretable ML for the grass-drought dataset
  - `Marker_Taxa_files/`: Folder containing plots to be displayed in `Marker_Taxa.md`
  - `Machine_Learning.ipynb`: Nested CV of Random Forest Classifier and interpretable ML with SHAP values with the grass-drought dataset
  - `data/`
    - `metadata.csv`: Table containing enrichment information per sample for the grass-drought dataset
    - `DADA2_ASV_count.Rdata`: Table containing ASV counts per sample of the grass-drought dataset
    - `DADA2_ASV_taxonomy.Rdata`: Table containing taxonomic annotation per ASV of the grass-drought dataset
    - `feature_tbl_{phylum/class/order/family/genus}.csv`: Table containing relative abundances of phyla/classes/orders/families/genera per sample with watering regime as target of the grass-drought dataset
  - `supplementary_analysis`:
    - `Logistic_Regression.py`: Python script running Nested CV of Logistic Regression Classifier with the grass-drought dataset
    - `Machine_Learning_Hold_Out.py`: Python script running Nested CV of Random Forest Classifier on the grass-drought dataset after the creation of a hold-out dataset from the grass-drought dataset and predicting on the hold-out dataset

The scripts need to be executed in the following order:

1. `Metagenomics_Analysis.Rmd`
2. `Machine_Learning.ipynb`
3. `Marker_Taxa.Rmd`

### Grass-Drought dataset:
**Naylor, D., DeGraaf, S., Purdom, E., Coleman-Derr, D.:**
Drought and host selection influence bacterial community dynamics in the grass root microbiome. The ISME Journal 11(12), 2691--2704 (2017)
<https://doi.org/10.1038/ismej.2017.118>
