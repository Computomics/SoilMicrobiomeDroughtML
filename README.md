# Interpretable Machine Learning Decodes Soil Microbiome’s Response to Drought Stress

This repository contains the code used for metagenomics analysis and ML applications for drought stress identification in the soil metagenome.

- `SoilMicrobiomeDroughtML/`
  - `Metagenomics_Analysis.Rmd`: Performs Data Processing, Diversity Analysis and Differential Abundance Analysis
  - `Marker_Taxa.Rmd`: Comparison of significant taxa from DAA and important taxa from interpretable ML
  - `Machine_Learning.ipynb`: Nested CV of Random Forest Classifier and interpretable ML with SHAP values
  - `data/`
    - `metadata.csv`: Table containing enrichment information per sample
    - `DADA2_ASV_count.Rdata`: Table containing ASV counts per sample
    - `DADA2_ASV_taxonomy.Rdata`: Table containing taxonomic annotation per ASV
    - `feature_tbl_phylum.csv`: Table containing relative abundances of phyla per sample with watering regime as target
    - `feature_tbl_class.csv`: Table containing relative abundances of classes per sample with watering regime as target
    - `feature_tbl_order.csv`: Table containing relative abundances of orders per sample with watering regime as target
    - `feature_tbl_family.csv`: Table containing relative abundances of families per sample with watering regime as target
    - `feature_tbl_genus.csv`: Table containing relative abundances of genera per sample with watering regime as target

The scripts need to be executed in the following order:

1. `Metagenomics_Analysis.Rmd`
2. `Machine_Learning.ipynb`
3. `Marker_Taxa.Rmd`

### Drought stress dataset:
**Naylor, D., DeGraaf, S., Purdom, E., Coleman-Derr, D.:**
Drought and host selection influence bacterial community dynamics in the grass root microbiome. The ISME Journal 11(12), 2691--2704 (2017)
<https://doi.org/10.1038/ismej.2017.118>
