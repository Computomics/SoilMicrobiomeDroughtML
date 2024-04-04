Metagenomics Analysis
================
Michelle Hagen

Perform 16S rRNA-based metagenomic analysis to the *‘drought-stress’ dataset* from
Naylor et al., 2017:

Naylor, D., DeGraaf, S., Purdom, E., Coleman-Derr, D.: Drought and host
selection influence bacterial community dynamics in the grass root
microbiome. The ISME Journal 11(12), 2691–2704 (2017)
<https://doi.org/10.1038/ismej.2017.118>

The script performs **data processing** steps and a **diversity
analysis** followed up with a **Differential Abundance Analysis** on ASV
level and all ranks (phylum to genus level).

``` r
library(phyloseq); print(paste0("phyloseq, version: ", packageVersion("phyloseq")))
```

    ## [1] "phyloseq, version: 1.44.0"

``` r
library(tidyverse); print(paste0("tidyverse, version: ", packageVersion("tidyverse")))
```

    ## [1] "tidyverse, version: 2.0.0"

``` r
library(ggpubr); print(paste0("ggpubr, version: ", packageVersion("ggpubr")))
```

    ## [1] "ggpubr, version: 0.6.0"

``` r
library(microbiomeMarker); print(paste0("microbiomeMarker, version: ", packageVersion("microbiomeMarker")))
```

    ## [1] "microbiomeMarker, version: 1.4.0"

``` r
library(ANCOMBC); print(paste0("ANCOMBC, version: ", packageVersion("ANCOMBC")))
```

    ## [1] "ANCOMBC, version: 2.0.2"

``` r
library(UpSetR); print(paste0("UpSetR, version: ", packageVersion("UpSetR")))
```

    ## [1] "UpSetR, version: 1.4.0"

**Load in Data:**

The **ASV tables** of the grass-drought dataset obtained from the DADA2
workflow are loaded in as RData objects directly into the workspace
(`data/DADA2_ASV_count.RData`, `data/DADA2_ASV_taxonomy.RData`) together
with the corresponding **metadata** file `data/metadata.csv`.

``` r
load("data/DADA2_ASV_count.RData")
load("data/DADA2_ASV_taxonomy.RData")
metadata <- read.table("data/metadata.csv", sep = ",", header = TRUE) %>% 
  column_to_rownames(var="Sample")
```

# Data Processing

## Filtering

ASVs that are not present in 95% or more of all samples are filtered
from the dataset. Samples with low readcounts (here, the 10% decile D1
of the dataset at 17,291 reads) are also discarded.

``` r
# Filter out ASVs that are not present in 95-100% of all samples and create new ASV_count and ASV_taxonomy tables
ASV_count_95 <- DADA2_ASV_count[,colSums(DADA2_ASV_count == 0) < (nrow(DADA2_ASV_count) * 0.95)]

t_ASV_taxonomy <- t(DADA2_ASV_taxonomy)
ASV_taxonomy_95 <- t(t_ASV_taxonomy[,colSums(DADA2_ASV_count == 0) < (nrow(DADA2_ASV_count) * 0.95)])
```

``` r
ps_filtered <- phyloseq(otu_table(ASV_count_95, taxa_are_rows=FALSE),
               sample_data(metadata),
               tax_table(ASV_taxonomy_95))

ps_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3276 taxa and 623 samples ]
    ## sample_data() Sample Data:       [ 623 samples by 3 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3276 taxa by 7 taxonomic ranks ]

``` r
D1 <- 17291

ps_filtered <- phyloseq::subset_samples(ps_filtered, phyloseq::sample_sums(ps_filtered) > D1)
(ps_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_filtered) > 0, ps_filtered))
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3276 taxa and 560 samples ]
    ## sample_data() Sample Data:       [ 560 samples by 3 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3276 taxa by 7 taxonomic ranks ]

## Rarefaction

The dataset is rarefied to the lowest sample size (here the 10 % decile
D1).

``` r
ps_rare_filtered <- phyloseq::rarefy_even_depth(ps_filtered, rngseed = 123, replace = FALSE)
```

    ## `set.seed(123)` was used to initialize repeatable random subsampling.

    ## Please record this for your records so others can reproduce.

    ## Try `set.seed(123); .Random.seed` for the full vector

    ## ...

``` r
ps_rare_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3276 taxa and 560 samples ]
    ## sample_data() Sample Data:       [ 560 samples by 3 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3276 taxa by 7 taxonomic ranks ]

Number of Control and Drought samples in the dataset:

``` r
table(sample_data(ps_rare_filtered)[["Watering_Regm"]])
```

    ## 
    ## Control Drought 
    ##     292     268

# Diversity Analysis

## Relative Abundances of Taxa

The top 10 most abundant taxa at each taxonomic rank (phylum - genus)
are displayed as stacked bar charts between the Control and Drought
groups to display visible differences in the relative abundances between
the two groups.

``` r
ps_ra_rare_filtered = phyloseq::transform_sample_counts(ps_rare_filtered, function(x){x / sum(x)})

ps_ra_top10_phy <- microbiomeutilities::aggregate_top_taxa2(ps_ra_rare_filtered, "Phylum", top = 10)
```

    ## Warning: replacing previous import 'ggplot2::alpha' by 'microbiome::alpha' when
    ## loading 'microbiomeutilities'

``` r
ps_ra_top10_cla <- microbiomeutilities::aggregate_top_taxa2(ps_ra_rare_filtered, "Class", top = 10)
ps_ra_top10_ord <- microbiomeutilities::aggregate_top_taxa2(ps_ra_rare_filtered, "Order", top = 10)
ps_ra_top10_fam <- microbiomeutilities::aggregate_top_taxa2(ps_ra_rare_filtered, "Family", top = 10)
ps_ra_top10_gen <- microbiomeutilities::aggregate_top_taxa2(ps_ra_rare_filtered, "Genus", top = 10)
ps_ra_top10_spe <- microbiomeutilities::aggregate_top_taxa2(ps_ra_rare_filtered, "Species", top = 10)
```

``` r
merge_and_transform <- function(ps_data, group_name) {
  ps_merge <- merge_samples(ps_data, group = group_name)
  sample_data(ps_merge)$Watering_Regm <- sample_names(ps_merge)
  ps_merge_comp <- transform_sample_counts(ps_merge, function(x) x / sum(x) * 100)
  return(ps_merge_comp)
}

ps_ra_top10_phy_merge <- merge_and_transform(ps_ra_top10_phy, "Watering_Regm")
ps_ra_top10_cla_merge <- merge_and_transform(ps_ra_top10_cla, "Watering_Regm")
ps_ra_top10_ord_merge <- merge_and_transform(ps_ra_top10_ord, "Watering_Regm")
ps_ra_top10_fam_merge <- merge_and_transform(ps_ra_top10_fam, "Watering_Regm")
ps_ra_top10_gen_merge <- merge_and_transform(ps_ra_top10_gen, "Watering_Regm")
```

``` r
top11_palette <- c("#290AD8", "#264DFF", "#3FA0FF", "#AAF7FF", "#B2FFB2", "#FFFFBF", "#FFE099", "#FFAD72", "#F76D5E", "#D82632", "#A50021")
```

``` r
# Phylum
p_ra_top10_Phylum <- plot_bar(ps_ra_top10_phy_merge, fill = "Phylum", x = "Watering_Regm") +
    geom_bar(aes(fill = Phylum),
             stat='identity',
             position='stack') +
  scale_fill_manual(values = top11_palette) +
  labs(x = "", y = "Relative Abundance\n", title = "Top10 Phyla") +
  theme_bw() +
  theme(legend.position="right") +
  theme(legend.title=element_blank())

p_ra_top10_Phylum$data$Phylum <- factor(p_ra_top10_Phylum$data$Phylum, levels = c("Acidobacteriota", "Actinobacteriota", "Bacteroidota", "Chloroflexi", "Firmicutes", "Gemmatimonadota", "Myxococcota", "Patescibacteria", "Proteobacteria", "Verrucomicrobiota", "Other"))

# Class
p_ra_top10_Class <- plot_bar(ps_ra_top10_cla_merge, fill = "Class", x = "Watering_Regm") +
  geom_bar(aes(fill = Class),  # use fill instead of color
           stat = 'identity',
           position = 'stack') +
  scale_fill_manual(values = top11_palette) +  # specify the new color palette
  labs(x = "", y = "Relative Abundance\n", title = "Top10 Classes") +
  theme_bw() +
  theme(legend.position = "right") +
  theme(legend.title = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

p_ra_top10_Class$data$Class <- factor(p_ra_top10_Class$data$Class, levels = c("Acidobacteriae", "Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidia", "Gammaproteobacteria", "Gemmatimonadetes", "Ktedonobacteria", "Polyangia", "Verrucomicrobiae", "Other"))

#Order
p_ra_top10_Order <- plot_bar(ps_ra_top10_ord_merge, fill = "Order", x = "Watering_Regm") +
    geom_bar(aes(fill = Order),
             stat='identity',
             position='stack') +
  scale_fill_manual(values = top11_palette) +
  labs(x = "", y = "Relative Abundance\n", title = "Top10 Orders") +
  theme_bw() +
  theme(legend.position="right") +
  theme(legend.title=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

p_ra_top10_Order$data$Order <- factor(p_ra_top10_Order$data$Order, levels = c("Acidobacteriales", "Bacillales", "Burkholderiales", "Chitinophagales", "Cytophagales", "Pedosphaerales", "Rhizobiales", "Sphingobacteriales", "Sphingomonadales", "Streptomycetales", "Other"))

#Family
p_ra_top10_Family <- plot_bar(ps_ra_top10_fam_merge, fill = "Family", x = "Watering_Regm") +
    geom_bar(aes(fill = Family),
             stat='identity',
             position='stack') +
  scale_fill_manual(values = top11_palette) +
  labs(x = "", y = "Relative Abundance\n", title = "Top10 Families") +
  theme_bw() +
  theme(legend.position="right") +
  theme(legend.title=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

p_ra_top10_Family$data$Family <- factor(p_ra_top10_Family$data$Family, levels = c("Acidobacteriaceae (Subgroup 1)", "Bacillaceae", "Chitinophagaceae", "Comamonadaceae", "Microscillaceae", "Pedosphaeraceae", "Sphingobacteriaceae", "Sphingomonadaceae", "Streptomycetaceae", "Other", "Unknown"))

#Genus
p_ra_top10_Genus <- plot_bar(ps_ra_top10_gen_merge, fill = "Genus", x = "Watering_Regm") +
    geom_bar(aes(fill = Genus),
             stat='identity',
             position='stack') +
  scale_fill_manual(values = top11_palette) +
  labs(x = "", y = "Relative Abundance\n", title = "Top10 Genera") +
  theme_bw() +
  theme(legend.position="right") +
  theme(legend.title=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

p_ra_top10_Genus$data$Genus <- factor(p_ra_top10_Genus$data$Genus, levels = c("Bacillus", "Bryobacter", "Candidatus Koribacter", "Candidatus Solibacter", "Gemmatimonas", "Mucilaginibacter", "Niastella", "Sphingomonas", "Streptomyces", "Other", "Unknown"))

# Arrange the plots using ggarrange
egg::ggarrange(p_ra_top10_Phylum, p_ra_top10_Class, p_ra_top10_Order, p_ra_top10_Family, p_ra_top10_Genus, ncol = 2)
```

<figure>
<img
src="Metagenomics_Analysis_files/figure-gfm/top10%20rel%20abundances-1.png"
alt="Relative Abundances per Rank of the Grass-Drought Dataset. Bar plots displaying the relative abundance of the top 10 taxa between the ’Control’ and ’Drought’ group on Phylum, Class, Order, Family and Genus level in alphabetical order." />
<figcaption aria-hidden="true"><strong>Relative Abundances per Rank of
the Grass-Drought Dataset.</strong> Bar plots displaying the relative
abundance of the top 10 taxa between the ’Control’ and ’Drought’ group
on Phylum, Class, Order, Family and Genus level in alphabetical
order.</figcaption>
</figure>

## Alpha Diversity

The variation within samples is displayed using the Shannon index
between the two watering regimes.

``` r
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(ps_rare_filtered, measures = "Shannon"),
  "Watering_Regm" = phyloseq::sample_data(ps_rare_filtered)$Watering_Regm
  )

wr <- list(c("Control", "Drought"))

adiv_wr <- adiv %>%
  gather(key = metric, value = value, "Shannon") %>%
  mutate(metric = factor(metric, levels = "Shannon")) %>%
  ggplot(aes(x = Watering_Regm, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Watering_Regm), height = 0, width = .2) +
  stat_compare_means(comparisons = wr, label = "p.signif")+
  labs(x = "", y = "") +
  scale_color_manual(values = c("#5a34e3", "#f54260")) +
  theme_minimal() +
  theme(legend.position="none") +
  #ggtitle("Alpha Diversity Shannon") +
  theme(plot.title = element_text(hjust = 0.5))

adiv_wr
```

<figure>
<img
src="Metagenomics_Analysis_files/figure-gfm/Alpha%20Diversity-1.png"
alt="Alpha Diversity for the Grass-Drought Dataset. Alpha diversity plot comparing the Control (blue) and ’Drought’ (red) watering regimes. (A) Boxplot of Shannon’s Diversity index for all samples comparing the ‘Control’ (blue) and ’Drought’ (red) watering regimes. Significance was determined using a non-parametric Wilcoxon rank sum test (* p &lt;0.05, ** p &lt;0.01, *** p &lt;0.001, **** p &lt;0.0001)." />
<figcaption aria-hidden="true"><strong>Alpha Diversity for the
Grass-Drought Dataset.</strong> Alpha diversity plot comparing the
Control (blue) and ’Drought’ (red) watering regimes. (A) Boxplot of
Shannon’s Diversity index for all samples comparing the ‘Control’ (blue)
and ’Drought’ (red) watering regimes. Significance was determined using
a non-parametric Wilcoxon rank sum test (* p &lt;0.05, ** p &lt;0.01,
*** p &lt;0.001, **** p &lt;0.0001).</figcaption>
</figure>

## Beta Diversity

The variation between samples is displayed using the Bray-Curtis
dissimilarities between the two watering regimes.

``` r
pcoa_bc = phyloseq::ordinate(ps_rare_filtered, "PCoA", "bray") 

pcoa_wr <- phyloseq::plot_ordination(ps_rare_filtered, pcoa_bc, color = "Watering_Regm") + 
  geom_point(size = 3) +
  #ggtitle("PCoA by Watering_Regm") +
  theme_minimal() +
  scale_color_manual(values = c("#5a34e3", "#f54260")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = 'bottom',
        legend.title = element_blank())

pcoa_wr
```

<figure>
<img src="Metagenomics_Analysis_files/figure-gfm/Beta%20diversity-1.png"
alt="Beta Diversity for the Grass-Drought Dataset. Principal Coordinate plot using Bray-Curtis dissimilarities colored by the ’Control’ (blue) and ’Drought’ (red) watering regimes." />
<figcaption aria-hidden="true"><strong>Beta Diversity for the
Grass-Drought Dataset.</strong> Principal Coordinate plot using
Bray-Curtis dissimilarities colored by the ’Control’ (blue) and
’Drought’ (red) watering regimes.</figcaption>
</figure>

## Differential Abundance Analysis

DAA was performed to find significant taxa for drought stress. Five
popular methods were tested out on ASV level and the three most suitable
methods were further used on all ranks (phylum-genus).

``` r
taxa_info <- data.frame(tax_table(ps_rare_filtered)) %>%
  rownames_to_column(var = "ASV")
```

## DAA on ASV Level

### Wilcoxon rank-sum test

``` r
ps_rare_filtered_clr <- microbiome::transform(ps_rare_filtered, "clr")

ps_wilcox_r <- data.frame(phyloseq::otu_table(ps_rare_filtered_clr))
ps_wilcox_r$Watering_Regm <- phyloseq::sample_data(ps_rare_filtered_clr)$Watering_Regm

wilcox_model <- function(df){
  wilcox.test(abund ~ Watering_Regm, data = df)
}
wilcox_pval <- function(df){
  wilcox.test(abund ~ Watering_Regm, data = df)$p.value
}


wilcox_results_r <- ps_wilcox_r %>%
  gather(key = ASV, value = abund, -Watering_Regm) %>%
  group_by(ASV) %>%
  nest() %>%
  mutate(wilcox_test = map(data, wilcox_model),
         p_value = map(data, wilcox_pval)) 

wilcox_results_r <- wilcox_results_r %>%
  dplyr::select(ASV, p_value) %>%
  unnest()
```

    ## Warning: `cols` is now required when using `unnest()`.
    ## ℹ Please use `cols = c(p_value)`.

``` r
sig_wilcox_r <- wilcox_results_r %>%
  full_join(taxa_info) %>%
  arrange(p_value) %>%
  mutate(BH_FDR = p.adjust(p_value, "BH")) %>%
  filter(BH_FDR < 0.05) %>%
  dplyr::select(ASV, p_value, BH_FDR, everything()) %>% 
  arrange(order(gtools::mixedorder(ASV)))
```

    ## Joining with `by = join_by(ASV)`

### edgeR

``` r
edger_microbiomeMarker <- run_edger(
  ps_rare_filtered,
  group = "Watering_Regm",
  method = "QLFT",
  taxa_rank = "none",
  transform = "identity",
  norm = "none",
  p_adjust = "BH",
  pvalue_cutoff = 0.05,
)

edger_marker <- marker_table(edger_microbiomeMarker) %>%
  as_tibble() %>% 
  arrange(order(gtools::mixedorder(feature))) %>% 
  dplyr::rename(ASV = feature) %>% 
  left_join(taxa_info, by = "ASV")
```

### DESeq2

``` r
deseq_microbiomeMarker <- run_deseq2(
  ps_rare_filtered,
  group = "Watering_Regm",
  confounders = character(0),
  contrast = NULL,
  taxa_rank = "none",  # To take only ASVs
  norm = "none",
  transform = "identity",
  fitType = "local",
  sfType = "poscounts",
  betaPrior = FALSE,
  useT = FALSE,
  p_adjust = "BH",
  pvalue_cutoff = 0.05
)
```

    ## converting counts to integer mode

``` r
deseq_marker <- marker_table(deseq_microbiomeMarker) %>%
  as_tibble() %>% 
  arrange(order(gtools::mixedorder(feature))) %>% 
  dplyr::rename(ASV = feature) %>% 
  left_join(taxa_info, by = "ASV")
```

### ALDEx2

``` r
aldex_microbiomeMarker <- run_aldex(
  ps_rare_filtered,
  group = "Watering_Regm",
  taxa_rank = "none",  # To take only ASVs
  transform = "identity",
  norm = "none",
  method = "wilcox.test",
  p_adjust = "BH",
  pvalue_cutoff = 0.05,
  mc_samples = 128,
  denom = "iqlr",
  paired = FALSE
)
```

    ## operating in serial mode

    ## computing iqlr centering

    ## New names:
    ## • `` -> `...1`
    ## • `` -> `...2`
    ## • `` -> `...3`
    ## • `` -> `...4`
    ## • `` -> `...5`
    ## • `` -> `...6`
    ## • `` -> `...7`
    ## • `` -> `...8`
    ## • `` -> `...9`
    ## • `` -> `...10`
    ## • `` -> `...11`
    ## • `` -> `...12`
    ## • `` -> `...13`
    ## • `` -> `...14`
    ## • `` -> `...15`
    ## • `` -> `...16`
    ## • `` -> `...17`
    ## • `` -> `...18`
    ## • `` -> `...19`
    ## • `` -> `...20`
    ## • `` -> `...21`
    ## • `` -> `...22`
    ## • `` -> `...23`
    ## • `` -> `...24`
    ## • `` -> `...25`
    ## • `` -> `...26`
    ## • `` -> `...27`
    ## • `` -> `...28`
    ## • `` -> `...29`
    ## • `` -> `...30`
    ## • `` -> `...31`
    ## • `` -> `...32`
    ## • `` -> `...33`
    ## • `` -> `...34`
    ## • `` -> `...35`
    ## • `` -> `...36`
    ## • `` -> `...37`
    ## • `` -> `...38`
    ## • `` -> `...39`
    ## • `` -> `...40`
    ## • `` -> `...41`
    ## • `` -> `...42`
    ## • `` -> `...43`
    ## • `` -> `...44`
    ## • `` -> `...45`
    ## • `` -> `...46`
    ## • `` -> `...47`
    ## • `` -> `...48`
    ## • `` -> `...49`
    ## • `` -> `...50`
    ## • `` -> `...51`
    ## • `` -> `...52`
    ## • `` -> `...53`
    ## • `` -> `...54`
    ## • `` -> `...55`
    ## • `` -> `...56`
    ## • `` -> `...57`
    ## • `` -> `...58`
    ## • `` -> `...59`
    ## • `` -> `...60`
    ## • `` -> `...61`
    ## • `` -> `...62`
    ## • `` -> `...63`
    ## • `` -> `...64`
    ## • `` -> `...65`
    ## • `` -> `...66`
    ## • `` -> `...67`
    ## • `` -> `...68`
    ## • `` -> `...69`
    ## • `` -> `...70`
    ## • `` -> `...71`
    ## • `` -> `...72`
    ## • `` -> `...73`
    ## • `` -> `...74`
    ## • `` -> `...75`
    ## • `` -> `...76`
    ## • `` -> `...77`
    ## • `` -> `...78`
    ## • `` -> `...79`
    ## • `` -> `...80`
    ## • `` -> `...81`
    ## • `` -> `...82`
    ## • `` -> `...83`
    ## • `` -> `...84`
    ## • `` -> `...85`
    ## • `` -> `...86`
    ## • `` -> `...87`
    ## • `` -> `...88`
    ## • `` -> `...89`
    ## • `` -> `...90`
    ## • `` -> `...91`
    ## • `` -> `...92`
    ## • `` -> `...93`
    ## • `` -> `...94`
    ## • `` -> `...95`
    ## • `` -> `...96`
    ## • `` -> `...97`
    ## • `` -> `...98`
    ## • `` -> `...99`
    ## • `` -> `...100`
    ## • `` -> `...101`
    ## • `` -> `...102`
    ## • `` -> `...103`
    ## • `` -> `...104`
    ## • `` -> `...105`
    ## • `` -> `...106`
    ## • `` -> `...107`
    ## • `` -> `...108`
    ## • `` -> `...109`
    ## • `` -> `...110`
    ## • `` -> `...111`
    ## • `` -> `...112`
    ## • `` -> `...113`
    ## • `` -> `...114`
    ## • `` -> `...115`
    ## • `` -> `...116`
    ## • `` -> `...117`
    ## • `` -> `...118`
    ## • `` -> `...119`
    ## • `` -> `...120`
    ## • `` -> `...121`
    ## • `` -> `...122`
    ## • `` -> `...123`
    ## • `` -> `...124`
    ## • `` -> `...125`
    ## • `` -> `...126`
    ## • `` -> `...127`
    ## • `` -> `...128`

``` r
aldex_marker <- marker_table(aldex_microbiomeMarker) %>%
  as_tibble() %>% 
  arrange(order(gtools::mixedorder(feature))) %>% 
  mutate(enrich_group = replace(enrich_group, enrich_group == "Control", "Watered")) %>%
  mutate(enrich_group = replace(enrich_group, enrich_group == "Drought", "Control")) %>% 
  mutate(enrich_group = replace(enrich_group, enrich_group == "Watered", "Drought")) %>% 
  dplyr::rename(ASV = feature) %>% 
  left_join(taxa_info, by = "ASV")
```

### ANCOM-BC2

``` r
ancom_da <- ancombc2(data = ps_rare_filtered, tax_level = NULL, fix_formula = "Watering_Regm",
              p_adj_method = "BH", lib_cut = 0, group = "Watering_Regm", struc_zero = FALSE,
              neg_lb = FALSE, alpha = 0.05, global = FALSE)
```

    ## Warning: The group variable has < 3 categories 
    ## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated

    ## Loading required package: foreach

    ## 
    ## Attaching package: 'foreach'

    ## The following objects are masked from 'package:purrr':
    ## 
    ##     accumulate, when

    ## Loading required package: rngtools

``` r
ancom_res <- data.frame(
  ASV = unlist(ancom_da$res$taxon),
  lfc = unlist(ancom_da$res$lfc_Watering_RegmDrought),
  se = unlist(ancom_da$res$lfc_Watering_RegmDrought),
  W = unlist(ancom_da$res$W_Watering_RegmDrought),
  p_val = unlist(ancom_da$res$p_Watering_RegmDrought),
  q_val = unlist(ancom_da$res$q_Watering_RegmDrought),
  diff = unlist(ancom_da$res$diff_Watering_RegmDrought))

ancombc2_marker <- ancom_res %>% 
  filter(q_val < 0.05) %>%  # adjusted p-values
  arrange(order(gtools::mixedorder(ASV))) %>%
  left_join(taxa_info) %>% 
  mutate(enrich_group = case_when(
    lfc >=0 ~ "Drought",
    lfc <0 ~ "Control"
  )) %>% 
  dplyr::rename(ef_ancombc2 = lfc) %>% 
  dplyr::rename(pvalue = p_val) %>% 
  dplyr::rename(padj = q_val) %>% 
  dplyr::select(ASV, enrich_group, ef_ancombc2, pvalue, padj) %>% 
  left_join(taxa_info, by = "ASV")
```

    ## Joining with `by = join_by(ASV)`

## DAA on all ranks

ANCOM-BC2, ALDEx2 and DESeq2 are selected as the top 3 best performing
tools for this dataset and therefore used on all taxonomic ranks
(phylum-genus).

``` r
Phylum_tax <- DADA2_ASV_taxonomy %>% 
  as_tibble() %>% 
  dplyr::select(Phylum) %>% 
  na.omit() %>%
  distinct()

Class_tax <- DADA2_ASV_taxonomy %>% 
  as_tibble() %>% 
  dplyr::select(Class) %>% 
  na.omit() %>%
  distinct()

Order_tax <- DADA2_ASV_taxonomy %>% 
  as_tibble() %>% 
  dplyr::select(Order) %>% 
  na.omit() %>%
  distinct()

Family_tax <- DADA2_ASV_taxonomy %>% 
  as_tibble() %>% 
  dplyr::select(Family) %>% 
  na.omit() %>%
  distinct()

Genus_tax <- DADA2_ASV_taxonomy %>% 
  as_tibble() %>% 
  dplyr::select(Genus) %>%
  na.omit() %>% 
  distinct()
```

### Phylum level DAA

``` r
# DESeq2
Phylum_deseq_microbiomeMarker <- run_deseq2(
  ps_rare_filtered,
  group = "Watering_Regm",
  confounders = character(0),
  contrast = NULL,
  taxa_rank = "Phylum",
  norm = "none",
  transform = "identity",
  fitType = "local",
  sfType = "poscounts",
  betaPrior = FALSE,
  useT = FALSE,
  p_adjust = "BH",
  pvalue_cutoff = 0.05
)
```

    ## converting counts to integer mode

``` r
Phylum_deseq_marker <- marker_table(Phylum_deseq_microbiomeMarker) %>%
  as_tibble() %>% 
  dplyr::filter(feature != "p__") %>% 
  arrange(order(gtools::mixedorder(enrich_group))) %>% 
  dplyr::rename(Phylum = feature) %>% 
  dplyr::rename(ef = ef_logFC) %>% 
  dplyr::rename(enrich = enrich_group) %>% 
  dplyr::rename(pvalue = pvalue) %>% 
  dplyr::rename(padj = padj) %>% 
  inner_join(Phylum_tax, by = "Phylum") %>% 
  mutate(tool="DESeq2")

Phylum_aldex_microbiomeMarker <- run_aldex(
  ps_rare_filtered,
  group = "Watering_Regm",
  taxa_rank = "Phylum",
  transform = "identity",
  norm = "none",
  method = "wilcox.test",
  p_adjust = "BH",
  pvalue_cutoff = 0.05,
  mc_samples = 128,
  denom = "iqlr",
  paired = FALSE
)
```

    ## operating in serial mode

    ## computing iqlr centering

    ## New names:
    ## • `` -> `...1`
    ## • `` -> `...2`
    ## • `` -> `...3`
    ## • `` -> `...4`
    ## • `` -> `...5`
    ## • `` -> `...6`
    ## • `` -> `...7`
    ## • `` -> `...8`
    ## • `` -> `...9`
    ## • `` -> `...10`
    ## • `` -> `...11`
    ## • `` -> `...12`
    ## • `` -> `...13`
    ## • `` -> `...14`
    ## • `` -> `...15`
    ## • `` -> `...16`
    ## • `` -> `...17`
    ## • `` -> `...18`
    ## • `` -> `...19`
    ## • `` -> `...20`
    ## • `` -> `...21`
    ## • `` -> `...22`
    ## • `` -> `...23`
    ## • `` -> `...24`
    ## • `` -> `...25`
    ## • `` -> `...26`
    ## • `` -> `...27`
    ## • `` -> `...28`
    ## • `` -> `...29`
    ## • `` -> `...30`
    ## • `` -> `...31`
    ## • `` -> `...32`
    ## • `` -> `...33`
    ## • `` -> `...34`
    ## • `` -> `...35`
    ## • `` -> `...36`
    ## • `` -> `...37`
    ## • `` -> `...38`
    ## • `` -> `...39`
    ## • `` -> `...40`
    ## • `` -> `...41`
    ## • `` -> `...42`
    ## • `` -> `...43`
    ## • `` -> `...44`
    ## • `` -> `...45`
    ## • `` -> `...46`
    ## • `` -> `...47`
    ## • `` -> `...48`
    ## • `` -> `...49`
    ## • `` -> `...50`
    ## • `` -> `...51`
    ## • `` -> `...52`
    ## • `` -> `...53`
    ## • `` -> `...54`
    ## • `` -> `...55`
    ## • `` -> `...56`
    ## • `` -> `...57`
    ## • `` -> `...58`
    ## • `` -> `...59`
    ## • `` -> `...60`
    ## • `` -> `...61`
    ## • `` -> `...62`
    ## • `` -> `...63`
    ## • `` -> `...64`
    ## • `` -> `...65`
    ## • `` -> `...66`
    ## • `` -> `...67`
    ## • `` -> `...68`
    ## • `` -> `...69`
    ## • `` -> `...70`
    ## • `` -> `...71`
    ## • `` -> `...72`
    ## • `` -> `...73`
    ## • `` -> `...74`
    ## • `` -> `...75`
    ## • `` -> `...76`
    ## • `` -> `...77`
    ## • `` -> `...78`
    ## • `` -> `...79`
    ## • `` -> `...80`
    ## • `` -> `...81`
    ## • `` -> `...82`
    ## • `` -> `...83`
    ## • `` -> `...84`
    ## • `` -> `...85`
    ## • `` -> `...86`
    ## • `` -> `...87`
    ## • `` -> `...88`
    ## • `` -> `...89`
    ## • `` -> `...90`
    ## • `` -> `...91`
    ## • `` -> `...92`
    ## • `` -> `...93`
    ## • `` -> `...94`
    ## • `` -> `...95`
    ## • `` -> `...96`
    ## • `` -> `...97`
    ## • `` -> `...98`
    ## • `` -> `...99`
    ## • `` -> `...100`
    ## • `` -> `...101`
    ## • `` -> `...102`
    ## • `` -> `...103`
    ## • `` -> `...104`
    ## • `` -> `...105`
    ## • `` -> `...106`
    ## • `` -> `...107`
    ## • `` -> `...108`
    ## • `` -> `...109`
    ## • `` -> `...110`
    ## • `` -> `...111`
    ## • `` -> `...112`
    ## • `` -> `...113`
    ## • `` -> `...114`
    ## • `` -> `...115`
    ## • `` -> `...116`
    ## • `` -> `...117`
    ## • `` -> `...118`
    ## • `` -> `...119`
    ## • `` -> `...120`
    ## • `` -> `...121`
    ## • `` -> `...122`
    ## • `` -> `...123`
    ## • `` -> `...124`
    ## • `` -> `...125`
    ## • `` -> `...126`
    ## • `` -> `...127`
    ## • `` -> `...128`

``` r
#ALDEx2
Phylum_aldex_marker <- marker_table(Phylum_aldex_microbiomeMarker) %>%
  as_tibble() %>% 
  arrange(order(gtools::mixedorder(enrich_group))) %>% 
  mutate(enrich_group = replace(enrich_group, enrich_group == "Control", "Watered")) %>%
  mutate(enrich_group = replace(enrich_group, enrich_group == "Drought", "Control")) %>% 
  mutate(enrich_group = replace(enrich_group, enrich_group == "Watered", "Drought")) %>% 
  dplyr::filter(feature != "p__") %>% 
  dplyr::rename(Phylum = feature) %>% 
  dplyr::rename(ef = ef_aldex) %>% 
  dplyr::rename(enrich = enrich_group) %>% 
  dplyr::rename(pvalue = pvalue) %>% 
  dplyr::rename(padj = padj)%>% 
  inner_join(Phylum_tax, by = "Phylum") %>% 
  mutate(tool="ALDEx2")

# ANCOM-BC2
Phylum_ancom_da <- ancombc2(
  data = ps_rare_filtered,
  tax_level = "Phylum",
  fix_formula = "Watering_Regm",
  p_adj_method = "BH",
  lib_cut = 0, group = "Watering_Regm",
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.05,
  global = FALSE)
```

    ## Warning: The group variable has < 3 categories 
    ## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated

``` r
Phylum_ancom_res <- data.frame(
  ASV = unlist(Phylum_ancom_da$res$taxon),
  lfc = unlist(Phylum_ancom_da$res$lfc_Watering_RegmDrought),  # log fold changes
  se = unlist(Phylum_ancom_da$res$lfc_Watering_RegmDrought),  # standard errors of lfc
  W = unlist(Phylum_ancom_da$res$W_Watering_RegmDrought),  # lfc/se
  p_val = unlist(Phylum_ancom_da$res$p_Watering_RegmDrought),  # p-values
  q_val = unlist(Phylum_ancom_da$res$q_Watering_RegmDrought),  # adjusted p-values
  diff = unlist(Phylum_ancom_da$res$diff_Watering_RegmDrought)) # TRUE if the taxon is significant (has q less than alpha 0.05)

Phylum_sig_ancombc2 <- Phylum_ancom_res %>% 
  filter(q_val < 0.05) %>%  # adjusted p-values
  arrange(order(gtools::mixedorder(ASV)))

Phylum_ancombc2_marker <- Phylum_sig_ancombc2 %>% 
  mutate(enrich_group = case_when(
    lfc >=0 ~ "Drought",
    lfc <0 ~ "Control"
  )) %>% 
  dplyr::rename(ef = lfc) %>% 
  dplyr::rename(pvalue = p_val) %>% 
  dplyr::rename(padj = q_val) %>% 
  mutate(across('ASV', str_replace, 'Phylum:', '')) %>% 
  dplyr::filter(ASV != "Kingdom:Bacteria") %>% 
  dplyr::rename(Phylum = ASV) %>% 
  dplyr::select(Phylum, enrich_group, ef, pvalue, padj) %>% 
  dplyr::rename(enrich = enrich_group) %>% 
  dplyr::rename(pvalue = pvalue) %>% 
  dplyr::rename(padj = padj) %>% 
  inner_join(Phylum_tax, by = "Phylum") %>% 
  mutate(tool="ANCOMBC2")
```

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `across("ASV", str_replace, "Phylum:", "")`.
    ## Caused by warning:
    ## ! The `...` argument of `across()` is deprecated as of dplyr 1.1.0.
    ## Supply arguments directly to `.fns` through an anonymous function instead.
    ## 
    ##   # Previously
    ##   across(a:b, mean, na.rm = TRUE)
    ## 
    ##   # Now
    ##   across(a:b, \(x) mean(x, na.rm = TRUE))

### Class level DAA

``` r
# DESeq2
Class_deseq_microbiomeMarker <- run_deseq2(
  ps_rare_filtered,
  group = "Watering_Regm",
  confounders = character(0),
  contrast = NULL,
  taxa_rank = "Class",
  norm = "none",
  transform = "identity",
  fitType = "local",
  sfType = "poscounts",
  betaPrior = FALSE,
  useT = FALSE,
  p_adjust = "BH",
  pvalue_cutoff = 0.05
)
```

    ## converting counts to integer mode

``` r
Class_deseq_marker <- marker_table(Class_deseq_microbiomeMarker) %>%
  as_tibble() %>% 
  arrange(order(gtools::mixedorder(enrich_group))) %>% 
  mutate(across('feature', str_replace, 'p__|_c__', '')) %>% 
  mutate(across('feature', str_replace, '^_c__$', 'test')) %>% 
  dplyr::filter(feature != "test") %>%
  dplyr::rename(Class = feature) %>% 
  dplyr::rename(ef = ef_logFC) %>% 
  dplyr::rename(enrich = enrich_group) %>% 
  dplyr::rename(pvalue = pvalue) %>% 
  dplyr::rename(padj = padj) %>% 
  inner_join(Class_tax, by = "Class") %>% 
  mutate(tool = "DESeq2")

# ALDEx2
Class_aldex_microbiomeMarker <- run_aldex(
  ps_rare_filtered,
  group = "Watering_Regm",
  taxa_rank = "Class",
  transform = "identity",
  norm = "none",
  method = "wilcox.test",
  p_adjust = "BH",
  pvalue_cutoff = 0.05,
  mc_samples = 128,
  denom = "iqlr",
  paired = FALSE
)
```

    ## operating in serial mode

    ## computing iqlr centering

    ## New names:
    ## • `` -> `...1`
    ## • `` -> `...2`
    ## • `` -> `...3`
    ## • `` -> `...4`
    ## • `` -> `...5`
    ## • `` -> `...6`
    ## • `` -> `...7`
    ## • `` -> `...8`
    ## • `` -> `...9`
    ## • `` -> `...10`
    ## • `` -> `...11`
    ## • `` -> `...12`
    ## • `` -> `...13`
    ## • `` -> `...14`
    ## • `` -> `...15`
    ## • `` -> `...16`
    ## • `` -> `...17`
    ## • `` -> `...18`
    ## • `` -> `...19`
    ## • `` -> `...20`
    ## • `` -> `...21`
    ## • `` -> `...22`
    ## • `` -> `...23`
    ## • `` -> `...24`
    ## • `` -> `...25`
    ## • `` -> `...26`
    ## • `` -> `...27`
    ## • `` -> `...28`
    ## • `` -> `...29`
    ## • `` -> `...30`
    ## • `` -> `...31`
    ## • `` -> `...32`
    ## • `` -> `...33`
    ## • `` -> `...34`
    ## • `` -> `...35`
    ## • `` -> `...36`
    ## • `` -> `...37`
    ## • `` -> `...38`
    ## • `` -> `...39`
    ## • `` -> `...40`
    ## • `` -> `...41`
    ## • `` -> `...42`
    ## • `` -> `...43`
    ## • `` -> `...44`
    ## • `` -> `...45`
    ## • `` -> `...46`
    ## • `` -> `...47`
    ## • `` -> `...48`
    ## • `` -> `...49`
    ## • `` -> `...50`
    ## • `` -> `...51`
    ## • `` -> `...52`
    ## • `` -> `...53`
    ## • `` -> `...54`
    ## • `` -> `...55`
    ## • `` -> `...56`
    ## • `` -> `...57`
    ## • `` -> `...58`
    ## • `` -> `...59`
    ## • `` -> `...60`
    ## • `` -> `...61`
    ## • `` -> `...62`
    ## • `` -> `...63`
    ## • `` -> `...64`
    ## • `` -> `...65`
    ## • `` -> `...66`
    ## • `` -> `...67`
    ## • `` -> `...68`
    ## • `` -> `...69`
    ## • `` -> `...70`
    ## • `` -> `...71`
    ## • `` -> `...72`
    ## • `` -> `...73`
    ## • `` -> `...74`
    ## • `` -> `...75`
    ## • `` -> `...76`
    ## • `` -> `...77`
    ## • `` -> `...78`
    ## • `` -> `...79`
    ## • `` -> `...80`
    ## • `` -> `...81`
    ## • `` -> `...82`
    ## • `` -> `...83`
    ## • `` -> `...84`
    ## • `` -> `...85`
    ## • `` -> `...86`
    ## • `` -> `...87`
    ## • `` -> `...88`
    ## • `` -> `...89`
    ## • `` -> `...90`
    ## • `` -> `...91`
    ## • `` -> `...92`
    ## • `` -> `...93`
    ## • `` -> `...94`
    ## • `` -> `...95`
    ## • `` -> `...96`
    ## • `` -> `...97`
    ## • `` -> `...98`
    ## • `` -> `...99`
    ## • `` -> `...100`
    ## • `` -> `...101`
    ## • `` -> `...102`
    ## • `` -> `...103`
    ## • `` -> `...104`
    ## • `` -> `...105`
    ## • `` -> `...106`
    ## • `` -> `...107`
    ## • `` -> `...108`
    ## • `` -> `...109`
    ## • `` -> `...110`
    ## • `` -> `...111`
    ## • `` -> `...112`
    ## • `` -> `...113`
    ## • `` -> `...114`
    ## • `` -> `...115`
    ## • `` -> `...116`
    ## • `` -> `...117`
    ## • `` -> `...118`
    ## • `` -> `...119`
    ## • `` -> `...120`
    ## • `` -> `...121`
    ## • `` -> `...122`
    ## • `` -> `...123`
    ## • `` -> `...124`
    ## • `` -> `...125`
    ## • `` -> `...126`
    ## • `` -> `...127`
    ## • `` -> `...128`

``` r
Class_aldex_marker <- marker_table(Class_aldex_microbiomeMarker) %>%
  as_tibble() %>% 
  arrange(order(gtools::mixedorder(enrich_group))) %>% 
  mutate(enrich_group = replace(enrich_group, enrich_group == "Control", "Watered")) %>%
  mutate(enrich_group = replace(enrich_group, enrich_group == "Drought", "Control")) %>% 
  mutate(enrich_group = replace(enrich_group, enrich_group == "Watered", "Drought")) %>% 
  dplyr::filter(feature != "c__") %>% 
  dplyr::rename(Class = feature) %>% 
  dplyr::rename(ef = ef_aldex) %>% 
  dplyr::rename(enrich = enrich_group) %>% 
  dplyr::rename(pvalue = pvalue) %>% 
  dplyr::rename(padj = padj) %>% 
  mutate(across('Class', str_replace, '_c__', '')) %>% 
  inner_join(Class_tax, by = "Class") %>% 
  mutate(tool = "ALDEx2")

# ANCOM-BC2
Class_ancom_da <- ancombc2(
  data = ps_rare_filtered,
  tax_level = "Class",
  fix_formula = "Watering_Regm",
  p_adj_method = "BH",
  lib_cut = 0, group = "Watering_Regm",
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.05,
  global = FALSE)
```

    ## Warning: The group variable has < 3 categories 
    ## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated

``` r
Class_ancom_res <- data.frame(
  ASV = unlist(Class_ancom_da$res$taxon),
  lfc = unlist(Class_ancom_da$res$lfc_Watering_RegmDrought),  # log fold changes
  se = unlist(Class_ancom_da$res$lfc_Watering_RegmDrought),  # standard errors of lfc
  W = unlist(Class_ancom_da$res$W_Watering_RegmDrought),  # lfc/se
  p_val = unlist(Class_ancom_da$res$p_Watering_RegmDrought),  # p-values
  q_val = unlist(Class_ancom_da$res$q_Watering_RegmDrought),  # adjusted p-values
  diff = unlist(Class_ancom_da$res$diff_Watering_RegmDrought)) # TRUE if the taxon is significant (has q less than alpha 0.05)

Class_sig_ancombc2 <- Class_ancom_res %>% 
  filter(q_val < 0.05) %>%  # adjusted p-values
  arrange(order(gtools::mixedorder(ASV)))

Class_ancombc2_marker <- Class_sig_ancombc2 %>% 
  mutate(enrich_group = case_when(
    lfc >=0 ~ "Drought",
    lfc <0 ~ "Control"
  )) %>% 
  dplyr::rename(ef = lfc) %>% 
  dplyr::rename(pvalue = p_val) %>% 
  dplyr::rename(padj = q_val) %>% 
  mutate(across('ASV', str_replace, 'Class:', '')) %>% 
  dplyr::filter(ASV != "Kingdom:Bacteria") %>% 
  dplyr::rename(Class = ASV) %>% 
  dplyr::select(Class, enrich_group, ef, pvalue, padj) %>% 
  dplyr::rename(enrich = enrich_group) %>% 
  dplyr::rename(pvalue = pvalue) %>% 
  dplyr::rename(padj = padj) %>% 
  mutate(across('Class', str_replace, 'Phylum:', '')) %>% 
  inner_join(Class_tax, by = "Class") %>% 
  mutate(tool = "ANCOMBC2")
```

### Order level DAA

``` r
#DESeq2
Order_deseq_microbiomeMarker <- run_deseq2(
  ps_rare_filtered,
  group = "Watering_Regm",
  confounders = character(0),
  contrast = NULL,
  taxa_rank = "Order",
  norm = "none",
  transform = "identity",
  fitType = "local",
  sfType = "poscounts",
  betaPrior = FALSE,
  useT = FALSE,
  p_adjust = "BH",
  pvalue_cutoff = 0.05
)
```

    ## converting counts to integer mode

``` r
Order_deseq_marker <- marker_table(Order_deseq_microbiomeMarker) %>%
  as_tibble() %>% 
  mutate(across('feature', str_replace, '_c___o__|_o__', '')) %>% 
  mutate(across('feature', str_replace, '^p__$', 'test')) %>% 
  dplyr::filter(feature != "test") %>% 
  arrange(order(gtools::mixedorder(enrich_group))) %>% 
  dplyr::rename(Order = feature) %>% 
  dplyr::rename(ef = ef_logFC) %>% 
  dplyr::rename(enrich = enrich_group) %>% 
  dplyr::rename(pvalue = pvalue) %>% 
  dplyr::rename(padj = padj) %>% 
  inner_join(Order_tax, by = "Order") %>% 
  mutate(tool = "DESeq2")

# ALDEx2
Order_aldex_microbiomeMarker <- run_aldex(
  ps_rare_filtered,
  group = "Watering_Regm",
  taxa_rank = "Order",
  transform = "identity",
  norm = "none",
  method = "wilcox.test",
  p_adjust = "BH",
  pvalue_cutoff = 0.05,
  mc_samples = 128,
  denom = "iqlr",
  paired = FALSE
)
```

    ## operating in serial mode

    ## computing iqlr centering

    ## New names:
    ## • `` -> `...1`
    ## • `` -> `...2`
    ## • `` -> `...3`
    ## • `` -> `...4`
    ## • `` -> `...5`
    ## • `` -> `...6`
    ## • `` -> `...7`
    ## • `` -> `...8`
    ## • `` -> `...9`
    ## • `` -> `...10`
    ## • `` -> `...11`
    ## • `` -> `...12`
    ## • `` -> `...13`
    ## • `` -> `...14`
    ## • `` -> `...15`
    ## • `` -> `...16`
    ## • `` -> `...17`
    ## • `` -> `...18`
    ## • `` -> `...19`
    ## • `` -> `...20`
    ## • `` -> `...21`
    ## • `` -> `...22`
    ## • `` -> `...23`
    ## • `` -> `...24`
    ## • `` -> `...25`
    ## • `` -> `...26`
    ## • `` -> `...27`
    ## • `` -> `...28`
    ## • `` -> `...29`
    ## • `` -> `...30`
    ## • `` -> `...31`
    ## • `` -> `...32`
    ## • `` -> `...33`
    ## • `` -> `...34`
    ## • `` -> `...35`
    ## • `` -> `...36`
    ## • `` -> `...37`
    ## • `` -> `...38`
    ## • `` -> `...39`
    ## • `` -> `...40`
    ## • `` -> `...41`
    ## • `` -> `...42`
    ## • `` -> `...43`
    ## • `` -> `...44`
    ## • `` -> `...45`
    ## • `` -> `...46`
    ## • `` -> `...47`
    ## • `` -> `...48`
    ## • `` -> `...49`
    ## • `` -> `...50`
    ## • `` -> `...51`
    ## • `` -> `...52`
    ## • `` -> `...53`
    ## • `` -> `...54`
    ## • `` -> `...55`
    ## • `` -> `...56`
    ## • `` -> `...57`
    ## • `` -> `...58`
    ## • `` -> `...59`
    ## • `` -> `...60`
    ## • `` -> `...61`
    ## • `` -> `...62`
    ## • `` -> `...63`
    ## • `` -> `...64`
    ## • `` -> `...65`
    ## • `` -> `...66`
    ## • `` -> `...67`
    ## • `` -> `...68`
    ## • `` -> `...69`
    ## • `` -> `...70`
    ## • `` -> `...71`
    ## • `` -> `...72`
    ## • `` -> `...73`
    ## • `` -> `...74`
    ## • `` -> `...75`
    ## • `` -> `...76`
    ## • `` -> `...77`
    ## • `` -> `...78`
    ## • `` -> `...79`
    ## • `` -> `...80`
    ## • `` -> `...81`
    ## • `` -> `...82`
    ## • `` -> `...83`
    ## • `` -> `...84`
    ## • `` -> `...85`
    ## • `` -> `...86`
    ## • `` -> `...87`
    ## • `` -> `...88`
    ## • `` -> `...89`
    ## • `` -> `...90`
    ## • `` -> `...91`
    ## • `` -> `...92`
    ## • `` -> `...93`
    ## • `` -> `...94`
    ## • `` -> `...95`
    ## • `` -> `...96`
    ## • `` -> `...97`
    ## • `` -> `...98`
    ## • `` -> `...99`
    ## • `` -> `...100`
    ## • `` -> `...101`
    ## • `` -> `...102`
    ## • `` -> `...103`
    ## • `` -> `...104`
    ## • `` -> `...105`
    ## • `` -> `...106`
    ## • `` -> `...107`
    ## • `` -> `...108`
    ## • `` -> `...109`
    ## • `` -> `...110`
    ## • `` -> `...111`
    ## • `` -> `...112`
    ## • `` -> `...113`
    ## • `` -> `...114`
    ## • `` -> `...115`
    ## • `` -> `...116`
    ## • `` -> `...117`
    ## • `` -> `...118`
    ## • `` -> `...119`
    ## • `` -> `...120`
    ## • `` -> `...121`
    ## • `` -> `...122`
    ## • `` -> `...123`
    ## • `` -> `...124`
    ## • `` -> `...125`
    ## • `` -> `...126`
    ## • `` -> `...127`
    ## • `` -> `...128`

``` r
Order_aldex_marker <- marker_table(Order_aldex_microbiomeMarker) %>%
  as_tibble() %>% 
  arrange(order(gtools::mixedorder(enrich_group))) %>% 
  mutate(enrich_group = replace(enrich_group, enrich_group == "Control", "Watered")) %>%
  mutate(enrich_group = replace(enrich_group, enrich_group == "Drought", "Control")) %>% 
  mutate(enrich_group = replace(enrich_group, enrich_group == "Watered", "Drought")) %>% 
  mutate(across('feature', str_replace, '_c___o__|_o__', '')) %>% 
  dplyr::rename(Order = feature) %>% 
  dplyr::rename(ef = ef_aldex) %>% 
  dplyr::rename(enrich = enrich_group) %>% 
  dplyr::rename(pvalue = pvalue) %>% 
  dplyr::rename(padj = padj) %>% 
  inner_join(Order_tax, by = "Order") %>% 
  mutate(tool = "ALDEx2")

# ANCOM-BC2
Order_ancom_da <- ancombc2(
  data = ps_rare_filtered,
  tax_level = "Order",
  fix_formula = "Watering_Regm",
  p_adj_method = "BH",
  lib_cut = 0, group = "Watering_Regm",
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.05,
  global = FALSE)
```

    ## Warning: The group variable has < 3 categories 
    ## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated

``` r
Order_ancom_res <- data.frame(
  ASV = unlist(Order_ancom_da$res$taxon),
  lfc = unlist(Order_ancom_da$res$lfc_Watering_RegmDrought),  # log fold changes
  se = unlist(Order_ancom_da$res$lfc_Watering_RegmDrought),  # standard errors of lfc
  W = unlist(Order_ancom_da$res$W_Watering_RegmDrought),  # lfc/se
  p_val = unlist(Order_ancom_da$res$p_Watering_RegmDrought),  # p-values
  q_val = unlist(Order_ancom_da$res$q_Watering_RegmDrought),  # adjusted p-values
  diff = unlist(Order_ancom_da$res$diff_Watering_RegmDrought)) # TRUE if the taxon is significant (has q less than alpha 0.05)

Order_sig_ancombc2 <- Order_ancom_res %>% 
  filter(q_val < 0.05) %>%  # adjusted p-values
  arrange(order(gtools::mixedorder(ASV)))

Order_ancombc2_marker <- Order_sig_ancombc2 %>% 
  mutate(enrich_group = case_when(
    lfc >=0 ~ "Drought",
    lfc <0 ~ "Control"
  )) %>% 
  dplyr::rename(ef = lfc) %>% 
  dplyr::rename(pvalue = p_val) %>% 
  dplyr::rename(padj = q_val) %>% 
  mutate(across('ASV', str_replace, 'Order:', '')) %>% 
  dplyr::filter(ASV != "Kingdom:Bacteria") %>% 
  dplyr::rename(Order = ASV) %>% 
  dplyr::select(Order, enrich_group, ef, pvalue, padj) %>% 
  dplyr::rename(enrich = enrich_group) %>% 
  dplyr::rename(pvalue = pvalue) %>% 
  dplyr::rename(padj = padj) %>% 
  mutate(across('Order', str_replace, 'Phylum:', '')) %>% 
  mutate(across('Order', str_replace, 'Class:', '')) %>% 
  inner_join(Order_tax, by = "Order") %>% 
  mutate(tool = "ANCOMBC2")
```

### Family level DAA

``` r
#DESeq2
Family_deseq_microbiomeMarker <- run_deseq2(
  ps_rare_filtered,
  group = "Watering_Regm",
  confounders = character(0),
  contrast = NULL,
  taxa_rank = "Family",
  norm = "none",
  transform = "identity",
  fitType = "local",
  sfType = "poscounts",
  betaPrior = FALSE,
  useT = FALSE,
  p_adjust = "BH",
  pvalue_cutoff = 0.05
)
```

    ## converting counts to integer mode

``` r
Family_deseq_marker <- marker_table(Family_deseq_microbiomeMarker) %>%
  as_tibble() %>% 
  mutate(across('feature', str_replace, '_c___o___f__|_o___f__|_f__', '')) %>% 
  mutate(across('feature', str_replace, '^p___c___o___f__$', 'test')) %>% 
  arrange(order(gtools::mixedorder(enrich_group))) %>% 
  mutate(across('feature', str_replace, '^p__$', 'test')) %>% 
  arrange(order(gtools::mixedorder(enrich_group))) %>% 
  dplyr::rename(Family = feature) %>% 
  dplyr::rename(ef = ef_logFC) %>% 
  dplyr::rename(enrich = enrich_group) %>% 
  dplyr::rename(pvalue = pvalue) %>% 
  dplyr::rename(padj = padj) %>% 
  inner_join(Family_tax, by = "Family") %>% 
  mutate(tool = "DESeq2")

# ALDEx2
Family_aldex_microbiomeMarker <- run_aldex(
  ps_rare_filtered,
  group = "Watering_Regm",
  taxa_rank = "Family",
  transform = "identity",
  norm = "none",
  method = "wilcox.test",
  p_adjust = "BH",
  pvalue_cutoff = 0.05,
  mc_samples = 128,
  denom = "iqlr",
  paired = FALSE
)
```

    ## operating in serial mode

    ## computing iqlr centering

    ## New names:
    ## • `` -> `...1`
    ## • `` -> `...2`
    ## • `` -> `...3`
    ## • `` -> `...4`
    ## • `` -> `...5`
    ## • `` -> `...6`
    ## • `` -> `...7`
    ## • `` -> `...8`
    ## • `` -> `...9`
    ## • `` -> `...10`
    ## • `` -> `...11`
    ## • `` -> `...12`
    ## • `` -> `...13`
    ## • `` -> `...14`
    ## • `` -> `...15`
    ## • `` -> `...16`
    ## • `` -> `...17`
    ## • `` -> `...18`
    ## • `` -> `...19`
    ## • `` -> `...20`
    ## • `` -> `...21`
    ## • `` -> `...22`
    ## • `` -> `...23`
    ## • `` -> `...24`
    ## • `` -> `...25`
    ## • `` -> `...26`
    ## • `` -> `...27`
    ## • `` -> `...28`
    ## • `` -> `...29`
    ## • `` -> `...30`
    ## • `` -> `...31`
    ## • `` -> `...32`
    ## • `` -> `...33`
    ## • `` -> `...34`
    ## • `` -> `...35`
    ## • `` -> `...36`
    ## • `` -> `...37`
    ## • `` -> `...38`
    ## • `` -> `...39`
    ## • `` -> `...40`
    ## • `` -> `...41`
    ## • `` -> `...42`
    ## • `` -> `...43`
    ## • `` -> `...44`
    ## • `` -> `...45`
    ## • `` -> `...46`
    ## • `` -> `...47`
    ## • `` -> `...48`
    ## • `` -> `...49`
    ## • `` -> `...50`
    ## • `` -> `...51`
    ## • `` -> `...52`
    ## • `` -> `...53`
    ## • `` -> `...54`
    ## • `` -> `...55`
    ## • `` -> `...56`
    ## • `` -> `...57`
    ## • `` -> `...58`
    ## • `` -> `...59`
    ## • `` -> `...60`
    ## • `` -> `...61`
    ## • `` -> `...62`
    ## • `` -> `...63`
    ## • `` -> `...64`
    ## • `` -> `...65`
    ## • `` -> `...66`
    ## • `` -> `...67`
    ## • `` -> `...68`
    ## • `` -> `...69`
    ## • `` -> `...70`
    ## • `` -> `...71`
    ## • `` -> `...72`
    ## • `` -> `...73`
    ## • `` -> `...74`
    ## • `` -> `...75`
    ## • `` -> `...76`
    ## • `` -> `...77`
    ## • `` -> `...78`
    ## • `` -> `...79`
    ## • `` -> `...80`
    ## • `` -> `...81`
    ## • `` -> `...82`
    ## • `` -> `...83`
    ## • `` -> `...84`
    ## • `` -> `...85`
    ## • `` -> `...86`
    ## • `` -> `...87`
    ## • `` -> `...88`
    ## • `` -> `...89`
    ## • `` -> `...90`
    ## • `` -> `...91`
    ## • `` -> `...92`
    ## • `` -> `...93`
    ## • `` -> `...94`
    ## • `` -> `...95`
    ## • `` -> `...96`
    ## • `` -> `...97`
    ## • `` -> `...98`
    ## • `` -> `...99`
    ## • `` -> `...100`
    ## • `` -> `...101`
    ## • `` -> `...102`
    ## • `` -> `...103`
    ## • `` -> `...104`
    ## • `` -> `...105`
    ## • `` -> `...106`
    ## • `` -> `...107`
    ## • `` -> `...108`
    ## • `` -> `...109`
    ## • `` -> `...110`
    ## • `` -> `...111`
    ## • `` -> `...112`
    ## • `` -> `...113`
    ## • `` -> `...114`
    ## • `` -> `...115`
    ## • `` -> `...116`
    ## • `` -> `...117`
    ## • `` -> `...118`
    ## • `` -> `...119`
    ## • `` -> `...120`
    ## • `` -> `...121`
    ## • `` -> `...122`
    ## • `` -> `...123`
    ## • `` -> `...124`
    ## • `` -> `...125`
    ## • `` -> `...126`
    ## • `` -> `...127`
    ## • `` -> `...128`

``` r
Family_aldex_marker <- marker_table(Family_aldex_microbiomeMarker) %>%
  as_tibble() %>% 
  arrange(order(gtools::mixedorder(enrich_group))) %>% 
  mutate(enrich_group = replace(enrich_group, enrich_group == "Control", "Watered")) %>%
  mutate(enrich_group = replace(enrich_group, enrich_group == "Drought", "Control")) %>% 
  mutate(enrich_group = replace(enrich_group, enrich_group == "Watered", "Drought")) %>% 
  mutate(across('feature', str_replace, '_o___f__|_f__|_c___o___f__', '')) %>%
  dplyr::rename(Family = feature) %>% 
  dplyr::rename(ef = ef_aldex) %>% 
  dplyr::rename(enrich = enrich_group) %>% 
  dplyr::rename(pvalue = pvalue) %>% 
  dplyr::rename(padj = padj) %>% 
  inner_join(Family_tax, by = "Family") %>% 
  mutate(tool = "ALDEx2")

# ANCOM-BC2
Family_ancom_da <- ancombc2(
  data = ps_rare_filtered,
  tax_level = "Family",
  fix_formula = "Watering_Regm",
  p_adj_method = "BH",
  lib_cut = 0, group = "Watering_Regm",
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.05,
  global = FALSE)
```

    ## Warning: The group variable has < 3 categories 
    ## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated

``` r
Family_ancom_res <- data.frame(
  ASV = unlist(Family_ancom_da$res$taxon),
  lfc = unlist(Family_ancom_da$res$lfc_Watering_RegmDrought),  # log fold changes
  se = unlist(Family_ancom_da$res$lfc_Watering_RegmDrought),  # standard errors of lfc
  W = unlist(Family_ancom_da$res$W_Watering_RegmDrought),  # lfc/se
  p_val = unlist(Family_ancom_da$res$p_Watering_RegmDrought),  # p-values
  q_val = unlist(Family_ancom_da$res$q_Watering_RegmDrought),  # adjusted p-values
  diff = unlist(Family_ancom_da$res$diff_Watering_RegmDrought)) # TRUE if the taxon is significant (has q less than alpha 0.05)

Family_sig_ancombc2 <- Family_ancom_res %>% 
  filter(q_val < 0.05) %>%  # adjusted p-values
  arrange(order(gtools::mixedorder(ASV)))

Family_ancombc2_marker <- Family_sig_ancombc2 %>% 
  mutate(enrich_group = case_when(
    lfc >=0 ~ "Drought",
    lfc <0 ~ "Control"
  )) %>% 
  dplyr::rename(ef = lfc) %>% 
  dplyr::rename(pvalue = p_val) %>% 
  dplyr::rename(padj = q_val) %>% 
  mutate(across('ASV', str_replace, 'Family:', '')) %>% 
  dplyr::filter(ASV != "Kingdom:Bacteria") %>% 
  dplyr::rename(Family = ASV) %>% 
  dplyr::select(Family, enrich_group, ef, pvalue, padj) %>% 
  dplyr::rename(enrich = enrich_group) %>% 
  dplyr::rename(pvalue = pvalue) %>% 
  dplyr::rename(padj = padj) %>%
  mutate(across('Family', str_replace, 'Phylum:', '')) %>% 
  mutate(across('Family', str_replace, 'Class:', '')) %>% 
  mutate(across('Family', str_replace, 'Order:', '')) %>% 
  inner_join(Family_tax, by = "Family") %>% 
  mutate(tool = "ANCOMBC2")
```

### Genus level DAA

``` r
#DESeq2
Genus_deseq_microbiomeMarker <- run_deseq2(
  ps_rare_filtered,
  group = "Watering_Regm",
  confounders = character(0),
  contrast = NULL,
  taxa_rank = "Genus",
  norm = "none",
  transform = "identity",
  fitType = "local",
  sfType = "poscounts",
  betaPrior = FALSE,
  useT = FALSE,
  p_adjust = "BH",
  pvalue_cutoff = 0.05
)
```

    ## converting counts to integer mode

``` r
Genus_deseq_marker <- marker_table(Genus_deseq_microbiomeMarker) %>%
  as_tibble() %>% 
  mutate(across('feature', str_replace, '_o___f___g__|_g__|_f___g__|_c__', '')) %>%
  mutate(across('feature', str_replace, '_o___f___g__', '')) %>%
  mutate(across('feature', str_replace, '^p___c___o___f___g__$', 'test')) %>% 
  mutate(across('feature', str_replace, '^p___c__$', 'test')) %>% 
  mutate(across('feature', str_replace, '^p___o___f___g__$', 'test')) %>% 
  dplyr::filter(feature != 'test') %>% 
  arrange(order(gtools::mixedorder(enrich_group))) %>% 
  dplyr::rename(Genus = feature) %>% 
  dplyr::rename(ef = ef_logFC) %>% 
  dplyr::rename(enrich = enrich_group) %>% 
  dplyr::rename(pvalue = pvalue) %>% 
  dplyr::rename(padj = padj) %>% 
  inner_join(Genus_tax, by = "Genus") %>% 
  mutate(tool = "DESeq2")

# ALDEx2
Genus_aldex_microbiomeMarker <- run_aldex(
  ps_rare_filtered,
  group = "Watering_Regm",
  taxa_rank = "Genus",
  transform = "identity",
  norm = "none",
  method = "wilcox.test",
  p_adjust = "BH",
  pvalue_cutoff = 0.05,
  mc_samples = 128,
  denom = "iqlr",
  paired = FALSE
)
```

    ## operating in serial mode

    ## computing iqlr centering

    ## New names:
    ## • `` -> `...1`
    ## • `` -> `...2`
    ## • `` -> `...3`
    ## • `` -> `...4`
    ## • `` -> `...5`
    ## • `` -> `...6`
    ## • `` -> `...7`
    ## • `` -> `...8`
    ## • `` -> `...9`
    ## • `` -> `...10`
    ## • `` -> `...11`
    ## • `` -> `...12`
    ## • `` -> `...13`
    ## • `` -> `...14`
    ## • `` -> `...15`
    ## • `` -> `...16`
    ## • `` -> `...17`
    ## • `` -> `...18`
    ## • `` -> `...19`
    ## • `` -> `...20`
    ## • `` -> `...21`
    ## • `` -> `...22`
    ## • `` -> `...23`
    ## • `` -> `...24`
    ## • `` -> `...25`
    ## • `` -> `...26`
    ## • `` -> `...27`
    ## • `` -> `...28`
    ## • `` -> `...29`
    ## • `` -> `...30`
    ## • `` -> `...31`
    ## • `` -> `...32`
    ## • `` -> `...33`
    ## • `` -> `...34`
    ## • `` -> `...35`
    ## • `` -> `...36`
    ## • `` -> `...37`
    ## • `` -> `...38`
    ## • `` -> `...39`
    ## • `` -> `...40`
    ## • `` -> `...41`
    ## • `` -> `...42`
    ## • `` -> `...43`
    ## • `` -> `...44`
    ## • `` -> `...45`
    ## • `` -> `...46`
    ## • `` -> `...47`
    ## • `` -> `...48`
    ## • `` -> `...49`
    ## • `` -> `...50`
    ## • `` -> `...51`
    ## • `` -> `...52`
    ## • `` -> `...53`
    ## • `` -> `...54`
    ## • `` -> `...55`
    ## • `` -> `...56`
    ## • `` -> `...57`
    ## • `` -> `...58`
    ## • `` -> `...59`
    ## • `` -> `...60`
    ## • `` -> `...61`
    ## • `` -> `...62`
    ## • `` -> `...63`
    ## • `` -> `...64`
    ## • `` -> `...65`
    ## • `` -> `...66`
    ## • `` -> `...67`
    ## • `` -> `...68`
    ## • `` -> `...69`
    ## • `` -> `...70`
    ## • `` -> `...71`
    ## • `` -> `...72`
    ## • `` -> `...73`
    ## • `` -> `...74`
    ## • `` -> `...75`
    ## • `` -> `...76`
    ## • `` -> `...77`
    ## • `` -> `...78`
    ## • `` -> `...79`
    ## • `` -> `...80`
    ## • `` -> `...81`
    ## • `` -> `...82`
    ## • `` -> `...83`
    ## • `` -> `...84`
    ## • `` -> `...85`
    ## • `` -> `...86`
    ## • `` -> `...87`
    ## • `` -> `...88`
    ## • `` -> `...89`
    ## • `` -> `...90`
    ## • `` -> `...91`
    ## • `` -> `...92`
    ## • `` -> `...93`
    ## • `` -> `...94`
    ## • `` -> `...95`
    ## • `` -> `...96`
    ## • `` -> `...97`
    ## • `` -> `...98`
    ## • `` -> `...99`
    ## • `` -> `...100`
    ## • `` -> `...101`
    ## • `` -> `...102`
    ## • `` -> `...103`
    ## • `` -> `...104`
    ## • `` -> `...105`
    ## • `` -> `...106`
    ## • `` -> `...107`
    ## • `` -> `...108`
    ## • `` -> `...109`
    ## • `` -> `...110`
    ## • `` -> `...111`
    ## • `` -> `...112`
    ## • `` -> `...113`
    ## • `` -> `...114`
    ## • `` -> `...115`
    ## • `` -> `...116`
    ## • `` -> `...117`
    ## • `` -> `...118`
    ## • `` -> `...119`
    ## • `` -> `...120`
    ## • `` -> `...121`
    ## • `` -> `...122`
    ## • `` -> `...123`
    ## • `` -> `...124`
    ## • `` -> `...125`
    ## • `` -> `...126`
    ## • `` -> `...127`
    ## • `` -> `...128`

``` r
Genus_aldex_marker <- marker_table(Genus_aldex_microbiomeMarker) %>%
  as_tibble() %>% 
  arrange(order(gtools::mixedorder(enrich_group))) %>% 
  mutate(enrich_group = replace(enrich_group, enrich_group == "Control", "Watered")) %>%
  mutate(enrich_group = replace(enrich_group, enrich_group == "Drought", "Control")) %>% 
  mutate(enrich_group = replace(enrich_group, enrich_group == "Watered", "Drought")) %>% 
  mutate(across('feature', str_replace, '^p___c___o___f___g__$', 'test')) %>% 
  dplyr::filter(feature != "test") %>% 
  mutate(across('feature', str_replace, '_c___o___f___g__|_o___f___g__|_f___g__|_g__', '')) %>%
  dplyr::rename(Genus = feature) %>% 
  dplyr::rename(ef = ef_aldex) %>% 
  dplyr::rename(enrich = enrich_group) %>% 
  dplyr::rename(pvalue = pvalue) %>% 
  dplyr::rename(padj = padj) %>% 
  inner_join(Genus_tax, by = "Genus") %>% 
  mutate(tool = "ALDEx2")

# ANCOM-BC2
Genus_ancom_da <- ancombc2(
  data = ps_rare_filtered,
  tax_level = "Genus",
  fix_formula = "Watering_Regm",
  p_adj_method = "BH",
  lib_cut = 0, group = "Watering_Regm",
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.05,
  global = FALSE)
```

    ## Warning: The group variable has < 3 categories 
    ## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated

``` r
Genus_ancom_res <- data.frame(
  ASV = unlist(Genus_ancom_da$res$taxon),
  lfc = unlist(Genus_ancom_da$res$lfc_Watering_RegmDrought),  # log fold changes
  se = unlist(Genus_ancom_da$res$lfc_Watering_RegmDrought),  # standard errors of lfc
  W = unlist(Genus_ancom_da$res$W_Watering_RegmDrought),  # lfc/se
  p_val = unlist(Genus_ancom_da$res$p_Watering_RegmDrought),  # p-values
  q_val = unlist(Genus_ancom_da$res$q_Watering_RegmDrought),  # adjusted p-values
  diff = unlist(Genus_ancom_da$res$diff_Watering_RegmDrought)) # TRUE if the taxon is significant (has q less than alpha 0.05)

Genus_sig_ancombc2 <- Genus_ancom_res %>% 
  filter(q_val < 0.05) %>%  # adjusted p-values
  arrange(order(gtools::mixedorder(ASV)))

Genus_ancombc2_marker <- Genus_sig_ancombc2 %>% 
  mutate(enrich_group = case_when(
    lfc >=0 ~ "Drought",
    lfc <0 ~ "Control"
  )) %>% 
  dplyr::rename(ef = lfc) %>% 
  dplyr::rename(pvalue = p_val) %>% 
  dplyr::rename(padj = q_val) %>% 
  mutate(across('ASV', str_replace, 'Genus:', '')) %>% 
  dplyr::filter(ASV != "Kingdom:Bacteria") %>% 
  dplyr::rename(Genus = ASV) %>% 
  dplyr::select(Genus, enrich_group, ef, pvalue, padj) %>% 
  dplyr::rename(enrich = enrich_group) %>% 
  dplyr::rename(pvalue = pvalue) %>% 
  dplyr::rename(padj = padj) %>%
  mutate(across('Genus', str_replace, 'Phylum:', '')) %>% 
  mutate(across('Genus', str_replace, 'Class:', '')) %>% 
  mutate(across('Genus', str_replace, 'Order:', '')) %>% 
  mutate(across('Genus', str_replace, 'Family:', '')) %>% 
  inner_join(Genus_tax, by = "Genus") %>% 
  mutate(tool = "ANCOMBC2")
```

``` r
Phylum_marker <- Phylum_deseq_marker %>% 
  bind_rows(Phylum_aldex_marker) %>% 
  bind_rows(Phylum_ancombc2_marker)

Class_marker <- Class_deseq_marker %>% 
  bind_rows(Class_aldex_marker) %>% 
  bind_rows(Class_ancombc2_marker)

Order_marker <- Order_deseq_marker %>% 
  bind_rows(Order_aldex_marker) %>% 
  bind_rows(Order_ancombc2_marker)

Family_marker <- Family_deseq_marker %>% 
  bind_rows(Family_aldex_marker) %>% 
  bind_rows(Family_ancombc2_marker)

Genus_marker <- Genus_deseq_marker %>% 
  bind_rows(Genus_aldex_marker) %>% 
  bind_rows(Genus_ancombc2_marker)
```

``` r
save(Phylum_marker, file = "data/Phylum_marker.RData")
save(Class_marker, file = "data/Class_marker.RData")
save(Order_marker, file = "data/Order_marker.RData")
save(Family_marker, file = "data/Family_marker.RData")
save(Genus_marker, file = "data/Genus_marker.RData")
```

## Display taxa/ASV intersection using UpSet plots

``` r
ASV_intersection <- list(Wilcox = sig_wilcox_r$ASV, DESeq2 = deseq_marker$ASV, edgeR = edger_marker$ASV, ALDEx2 = aldex_marker$ASV, ANCOMBC2 = ancombc2_marker$ASV)
top3_ASV_intersection <- list(DESeq2 = deseq_marker$ASV, ALDEx2 = aldex_marker$ASV, ANCOMBC2 = ancombc2_marker$ASV)
Phylum_intersection <- list(DESeq2 = Phylum_deseq_marker$Phylum, ALDEx2 = Phylum_aldex_marker$Phylum, ANCOMBC2 = Phylum_ancombc2_marker$Phylum)
Class_intersection <- list(DESeq2 = Class_deseq_marker$Class, ALDEx2 = Class_aldex_marker$Class, ANCOMBC2 = Class_ancombc2_marker$Class)
Order_intersection <- list(DESeq2 = Order_deseq_marker$Order, ALDEx2 = Order_aldex_marker$Order, ANCOMBC2 = Order_ancombc2_marker$Order)
Family_intersection <- list(DESeq2 = Family_deseq_marker$Family, ALDEx2 = Family_aldex_marker$Family, ANCOMBC2 = Family_ancombc2_marker$Family)
Genus_intersection <- list(DESeq2 = Genus_deseq_marker$Genus, ALDEx2 = Genus_aldex_marker$Genus, ANCOMBC2 = Genus_ancombc2_marker$Genus)
```

``` r
aldex_col <- "#bb0a21"
ancombc_col <- "#280659"
deseq_col <- "#faaf40"
edger_col <- "#014700"
wilcox_col <- "#d30085"

drought_col <- "#f54260"
control_col <- "#5a34e3"

soil_col <- "#ffd300"
root_col <- "#a35c00"
rhizosphere_col <- "#007355"
```

``` r
upset(fromList(ASV_intersection), order.by = "freq", point.size=4, sets.bar.color=c(edger_col, wilcox_col, ancombc_col, deseq_col, aldex_col), mainbar.y.label = "Tool Intersections", sets.x.label = "ASVs per Tool", text.scale = c(1.5, 1, 1, 1, 1.5, 1.2))
```

<figure>
<img
src="Metagenomics_Analysis_files/figure-gfm/UpsetPlot%20between%20different%20DAA%20Methods-1.png"
alt="ASV Intersections between different DAA Tools on ASV level of the Grass-Drought Dataset. Upset plots displaying the overlap and uniqueness of significant taxa on ASV level identified by DAA methods (ALDEx2, DESeq2, ANCOM-BC2, non-parametric Wilcoxon rank-sum test, and edgeR). The horizontal bars show the total number of taxa for each tool, while the vertical bars show the number of shared taxa between corresponding sets, sorted by the total number of shared taxa. All tools use an alpha threshold of 0.05 for significance" />
<figcaption aria-hidden="true"><strong>ASV Intersections between
different DAA Tools on ASV level of the Grass-Drought Dataset.</strong>
Upset plots displaying the overlap and uniqueness of significant taxa on
ASV level identified by DAA methods (ALDEx2, DESeq2, ANCOM-BC2,
non-parametric Wilcoxon rank-sum test, and edgeR). The horizontal bars
show the total number of taxa for each tool, while the vertical bars
show the number of shared taxa between corresponding sets, sorted by the
total number of shared taxa. All tools use an alpha threshold of 0.05
for significance</figcaption>
</figure>

``` r
upset(fromList(Phylum_intersection), order.by = "freq", point.size=4, sets.bar.color=c(ancombc_col, deseq_col, aldex_col), mainbar.y.label = "Tool Intersections", sets.x.label = "Phyla per Tool", text.scale = c(1.5, 1, 1, 1, 1.5, 1.2))
```

<figure>
<img
src="Metagenomics_Analysis_files/figure-gfm/Phylum%20UpsetPlot%20between%20top3%20DAA%20Methods-1.png"
alt="Phylum Intersections between different DAA Tools on the Grass-Drought Dataset. Upset plots displaying the overlap and uniqueness of significant phyla identified by top 3 DAA methods (ALDEx2, DESeq2, ANCOM-BC2). The horizontal bars show the total number of phyla for each tool, while the vertical bars show the number of shared phyla between corresponding sets, sorted by the total number of shared phyla. All tools use an alpha threshold of 0.05 for significance" />
<figcaption aria-hidden="true"><strong>Phylum Intersections between
different DAA Tools on the Grass-Drought Dataset.</strong> Upset plots
displaying the overlap and uniqueness of significant phyla identified by
top 3 DAA methods (ALDEx2, DESeq2, ANCOM-BC2). The horizontal bars show
the total number of phyla for each tool, while the vertical bars show
the number of shared phyla between corresponding sets, sorted by the
total number of shared phyla. All tools use an alpha threshold of 0.05
for significance</figcaption>
</figure>

``` r
upset(fromList(Class_intersection), order.by = "freq", point.size=4, sets.bar.color=c(ancombc_col, deseq_col, aldex_col), mainbar.y.label = "Tool Intersections", sets.x.label = "Classes per Tool", text.scale = c(1.5, 1, 1, 1, 1.5, 1.2))
```

<figure>
<img
src="Metagenomics_Analysis_files/figure-gfm/Class%20UpsetPlot%20between%20top3%20DAA%20Methods-1.png"
alt="Class Intersections between different DAA Tools on the Grass-Drought Dataset. Upset plots displaying the overlap and uniqueness of significant classes identified by top 3 DAA methods (ALDEx2, DESeq2, ANCOM-BC2). The horizontal bars show the total number of classes for each tool, while the vertical bars show the number of shared classes between corresponding sets, sorted by the total number of shared classes. All tools use an alpha threshold of 0.05 for significance" />
<figcaption aria-hidden="true"><strong>Class Intersections between
different DAA Tools on the Grass-Drought Dataset.</strong> Upset plots
displaying the overlap and uniqueness of significant classes identified
by top 3 DAA methods (ALDEx2, DESeq2, ANCOM-BC2). The horizontal bars
show the total number of classes for each tool, while the vertical bars
show the number of shared classes between corresponding sets, sorted by
the total number of shared classes. All tools use an alpha threshold of
0.05 for significance</figcaption>
</figure>

``` r
upset(fromList(Order_intersection), order.by = "freq", point.size=4, sets.bar.color=c(ancombc_col, deseq_col, aldex_col), mainbar.y.label = "Tool Intersections", sets.x.label = "Orders per Tool", text.scale = c(1.5, 1, 1, 1, 1.5, 1.2))
```

<figure>
<img
src="Metagenomics_Analysis_files/figure-gfm/Order%20UpsetPlot%20between%20top3%20DAA%20Methods-1.png"
alt="Order Intersections between different DAA Tools on the Grass-Drought Dataset. Upset plots displaying the overlap and uniqueness of significant orders identified by top 3 DAA methods (ALDEx2, DESeq2, ANCOM-BC2). The horizontal bars show the total number of orders for each tool, while the vertical bars show the number of shared orders between corresponding sets, sorted by the total number of shared orders. All tools use an alpha threshold of 0.05 for significance" />
<figcaption aria-hidden="true"><strong>Order Intersections between
different DAA Tools on the Grass-Drought Dataset.</strong> Upset plots
displaying the overlap and uniqueness of significant orders identified
by top 3 DAA methods (ALDEx2, DESeq2, ANCOM-BC2). The horizontal bars
show the total number of orders for each tool, while the vertical bars
show the number of shared orders between corresponding sets, sorted by
the total number of shared orders. All tools use an alpha threshold of
0.05 for significance</figcaption>
</figure>

``` r
upset(fromList(Family_intersection), order.by = "freq", point.size=4, sets.bar.color=c(ancombc_col, deseq_col, aldex_col), mainbar.y.label = "Tool Intersections", sets.x.label = "Families per Tool", text.scale = c(1.5, 1, 1, 1, 1.5, 1.2))
```

<figure>
<img
src="Metagenomics_Analysis_files/figure-gfm/Family%20UpsetPlot%20between%20top3%20DAA%20Methods-1.png"
alt="Family Intersections between different DAA Tools on the Grass-Drought Dataset. Upset plots displaying the overlap and uniqueness of significant families identified by top 3 DAA methods (ALDEx2, DESeq2, ANCOM-BC2). The horizontal bars show the total number of families for each tool, while the vertical bars show the number of shared families between corresponding sets, sorted by the total number of shared families. All tools use an alpha threshold of 0.05 for significance" />
<figcaption aria-hidden="true"><strong>Family Intersections between
different DAA Tools on the Grass-Drought Dataset.</strong> Upset plots
displaying the overlap and uniqueness of significant families identified
by top 3 DAA methods (ALDEx2, DESeq2, ANCOM-BC2). The horizontal bars
show the total number of families for each tool, while the vertical bars
show the number of shared families between corresponding sets, sorted by
the total number of shared families. All tools use an alpha threshold of
0.05 for significance</figcaption>
</figure>

``` r
upset(fromList(Genus_intersection), order.by = "freq", point.size=4, sets.bar.color=c(ancombc_col, deseq_col, aldex_col), mainbar.y.label = "Tool Intersections", sets.x.label = "Genera per Tool", text.scale = c(1.5, 1, 1, 1, 1.5, 1.2))
```

<figure>
<img
src="Metagenomics_Analysis_files/figure-gfm/Genus%20UpsetPlot%20between%20top3%20DAA%20Methods-1.png"
alt="Genus Intersections between different DAA Tools on the Grass-Drought Dataset. Upset plots displaying the overlap and uniqueness of significant genera identified by top 3 DAA methods (ALDEx2, DESeq2, ANCOM-BC2). The horizontal bars show the total number of genera for each tool, while the vertical bars show the number of shared genera between corresponding sets, sorted by the total number of shared genera. All tools use an alpha threshold of 0.05 for significance" />
<figcaption aria-hidden="true"><strong>Genus Intersections between
different DAA Tools on the Grass-Drought Dataset.</strong> Upset plots
displaying the overlap and uniqueness of significant genera identified
by top 3 DAA methods (ALDEx2, DESeq2, ANCOM-BC2). The horizontal bars
show the total number of genera for each tool, while the vertical bars
show the number of shared genera between corresponding sets, sorted by
the total number of shared genera. All tools use an alpha threshold of
0.05 for significance</figcaption>
</figure>
