Marker Taxa
================
Michelle Hagen

Creating plots displaying **enrichment of taxa** towards **Drought** or
**Control** conditions for ANCOM-BC2, ALDEx2, DESeq2 and the SHAP value
enrichment results combined with the corresponding **significance** or
**importance** per taxon for each rank (**phylum** to **genus**).

``` r
library(tidyverse); print(paste0("tidyverse, version: ", packageVersion("tidyverse")))
```

    ## [1] "tidyverse, version: 2.0.0"

``` r
library(phyloseq); print(paste0("phyloseq, version: ", packageVersion("phyloseq")))
```

    ## [1] "phyloseq, version: 1.44.0"

``` r
library(egg); print(paste0("egg, version: ", packageVersion("egg")))
```

    ## [1] "egg, version: 0.4.5"

## Calculate SHAP enrichment

The contribution of each taxon towards prediction of drought stress has
been calculated for each fold of the cv in the `Machine_Learning.ipynb`
script. Only if there is consistency in the enrichment pattern (4/5 or
5/5 consistency), the enrichment trend is used for the taxon.

``` r
# enrichment information of 4/5 or 5/5 folds needs to be matching for a final trend
threshold <-  4
```

``` r
load("data/Phylum_marker.RData")
load("data/Class_marker.RData")
load("data/Order_marker.RData")
load("data/Family_marker.RData")
load("data/Genus_marker.RData")

all_Phylum <- Phylum_marker %>% dplyr::select(Phylum) %>% distinct()
all_Class <- Class_marker %>% dplyr::select(Class) %>% distinct()
all_Order <- Order_marker %>% dplyr::select(Order) %>% distinct()
all_Family <- Family_marker %>% dplyr::select(Family) %>% distinct()
all_Genus <- Genus_marker %>% dplyr::select(Genus) %>% distinct()

Genus_marker <- Genus_marker %>% mutate(Genus = ifelse(Genus == "Burkholderia-Caballeronia-Paraburkholderia", "Burkholderia", Genus))
```

``` r
Phylum_SHAP <- read.table("data/SHAP_feature_importance_phylum.csv", sep = ",", header = TRUE) %>% 
  dplyr::rename(Taxon = features)

Class_SHAP <- read.table("data/SHAP_feature_importance_class.csv", sep = ",", header = TRUE) %>% 
  dplyr::rename(Taxon = features)

Order_SHAP <- read.table("data/SHAP_feature_importance_order.csv", sep = ",", header = TRUE) %>% 
  dplyr::rename(Taxon = features)

Family_SHAP <- read.table("data/SHAP_feature_importance_family.csv", sep = ",", header = TRUE) %>% 
  dplyr::rename(Taxon = features)

Genus_SHAP <- read.table("data/SHAP_feature_importance_genus.csv", sep = ",", header = TRUE) %>% 
  dplyr::rename(Taxon = features) %>% 
  mutate(Taxon = ifelse(Taxon == "Burkholderia-Caballeronia-Paraburkholderia", "Burkholderia", Taxon))
```

``` r
read_and_process_csv <- function(file_path, join_data, name_column) {
  df <- read.table(file_path, sep = ",", header = TRUE) %>%
    column_to_rownames("Taxa") %>%
    mutate(countDrought = rowSums(. == "Drought")) %>%
    mutate(Consensus = case_when(
      ((threshold == 3) & (countDrought %in% 3:5)) ~ "Drought",
      ((threshold == 3) & (countDrought %in% 0:2)) ~ "Control",
      ((threshold == 4) & (countDrought %in% 4:5)) ~ "Drought",
      ((threshold == 4) & (countDrought %in% 0:1)) ~ "Control"
    )) %>%
    dplyr::select(Consensus) %>%
    dplyr::rename(enrich = Consensus) %>%
    mutate(tool = "SHAP") %>%
    rownames_to_column(name_column) %>%
    inner_join(join_data)
  
  return(df)
}

Phylum_SHAP_enrich <- read_and_process_csv("data/SHAP_enriched_phylum.csv", all_Phylum, "Phylum") %>% 
  dplyr::rename(Taxon = Phylum)
```

    ## Joining with `by = join_by(Phylum)`

``` r
Class_SHAP_enrich <- read_and_process_csv("data/SHAP_enriched_class.csv", all_Class, "Class") %>% 
  dplyr::rename(Taxon = Class)
```

    ## Joining with `by = join_by(Class)`

``` r
Order_SHAP_enrich <- read_and_process_csv("data/SHAP_enriched_order.csv", all_Order, "Order") %>% 
  dplyr::rename(Taxon = Order)
```

    ## Joining with `by = join_by(Order)`

``` r
Family_SHAP_enrich <- read_and_process_csv("data/SHAP_enriched_family.csv", all_Family, "Family") %>% 
  dplyr::rename(Taxon = Family)
```

    ## Joining with `by = join_by(Family)`

``` r
Genus_SHAP_enrich <- read_and_process_csv("data/SHAP_enriched_genus.csv", all_Genus, "Genus") %>%
  mutate(Genus = ifelse(Genus == "Burkholderia-Caballeronia-Paraburkholderia", "Burkholderia", Genus)) %>% 
  dplyr::rename(Taxon = Genus)
```

    ## Joining with `by = join_by(Genus)`

The top25 significant taxa of ANCOM-BC2 are selected for plotting. (Top
15 on Phylum rank.)

``` r
ANCOMBC2_top25 <- function(data, tool_name, taxon_name) {
  result <- data %>%
    dplyr::filter(tool == tool_name) %>%
    dplyr::rename(Taxon = {{taxon_name}}) %>%
    mutate(log_padj = -log10(padj)) %>%
    dplyr::select(Taxon, log_padj) %>%
    dplyr::top_n(25, wt = log_padj) %>% 
    dplyr::arrange(desc(log_padj))
  return(result)
}

ANCOMBC2_top25_Phylum <- ANCOMBC2_top25(Phylum_marker, "ANCOMBC2", "Phylum")
ANCOMBC2_top25_Class <- ANCOMBC2_top25(Class_marker, "ANCOMBC2", "Class")
ANCOMBC2_top25_Order <- ANCOMBC2_top25(Order_marker, "ANCOMBC2", "Order")
ANCOMBC2_top25_Family <- ANCOMBC2_top25(Family_marker, "ANCOMBC2", "Family")
ANCOMBC2_top25_Genus <- ANCOMBC2_top25(Genus_marker, "ANCOMBC2", "Genus") %>%
  mutate(Taxon = ifelse(Taxon == "Burkholderia-Caballeronia-Paraburkholderia", "Burkholderia", Taxon))
```

``` r
Phylum_data <- Phylum_marker %>% 
  dplyr::rename(Taxon = Phylum) %>% 
  bind_rows(Phylum_SHAP_enrich) %>% 
  inner_join(ANCOMBC2_top25_Phylum, by = "Taxon") %>% 
  inner_join(Phylum_SHAP, by = "Taxon") %>% 
  arrange(desc(log_padj)) %>% 
  drop_na(enrich)

Class_data <- Class_marker %>% 
  dplyr::rename(Taxon = Class) %>% 
  na.omit() %>% 
  bind_rows(Class_SHAP_enrich) %>% 
  inner_join(ANCOMBC2_top25_Class, by = "Taxon") %>% 
  inner_join(Class_SHAP, by = "Taxon") %>% 
  arrange(desc(log_padj)) %>% 
  drop_na(enrich)

Order_data <- Order_marker %>% 
  dplyr::rename(Taxon = Order) %>% 
  na.omit() %>% 
  bind_rows(Order_SHAP_enrich) %>% 
  inner_join(ANCOMBC2_top25_Order, by = "Taxon") %>% 
  inner_join(Order_SHAP, by = "Taxon") %>% 
  arrange(desc(log_padj)) %>% 
  drop_na(enrich)

Family_data <- Family_marker %>% 
  dplyr::rename(Taxon = Family) %>% 
  na.omit() %>% 
  bind_rows(Family_SHAP_enrich) %>% 
  inner_join(ANCOMBC2_top25_Family, by = "Taxon") %>% 
  inner_join(Family_SHAP, by = "Taxon") %>% 
  arrange(desc(log_padj)) %>% 
  drop_na(enrich)

Genus_data <- Genus_marker %>% 
  dplyr::rename(Taxon = Genus) %>% 
  na.omit() %>% 
  bind_rows(Genus_SHAP_enrich) %>% 
  inner_join(ANCOMBC2_top25_Genus, by = "Taxon") %>% 
  inner_join(Genus_SHAP, by = "Taxon") %>% 
  arrange(desc(log_padj)) %>% 
  drop_na(enrich)
```

## Create Heatmaps + Bar plots for Taxon Enrichment/Significance/Importance

``` r
generate_heatmap_and_barplots <- function(data) {
  data <- data %>% rename(Enrichment = enrich)
  
  custom_order <- data %>%
    select(Taxon, log_padj) %>%
    distinct() %>%
    pull(Taxon)

  tool_enrich <- data %>%
    mutate(Taxon = factor(Taxon, levels = rev(custom_order)))

  sign_imp <- data %>%
    select(Taxon, log_padj, mean_shap) %>%
    distinct()

  hm <- ggplot(data = tool_enrich, aes(x = factor(tool, levels = c("ANCOMBC2", "ALDEx2", "DESeq2", "SHAP")), y = Taxon, fill = Enrichment)) +
    geom_tile(height = 0.9, width = 0.9) +
    coord_fixed() +
    theme_minimal(base_size = 10) +
    scale_x_discrete(position = "bottom", guide = guide_axis(angle = 45)) +
    scale_fill_manual(values = c("#5a34e3", "#f54260"), na.value = "white") +
    theme(legend.position = 'bottom',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90),
          axis.line.x = element_line(color = 'black'),
          axis.ticks.x = element_line(color = 'black'))

  bp1 <- ggplot(data = sign_imp, aes(x = log_padj, y = reorder(Taxon, log_padj))) +
    geom_bar(stat = 'identity', width = 0.9, fill = 'orange') +
    theme_minimal(base_size = 10) +
    labs(x = expression(-log[10]("p_adjust-value"))) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 0),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          axis.line.x = element_line(color = 'black'),
          axis.ticks.x = element_line(color = 'black'))

  bp2 <- ggplot(data = sign_imp, aes(x = mean_shap, y = reorder(Taxon, log_padj))) +
    geom_bar(stat = 'identity', width = 0.9, fill = 'darkgreen') +
    theme_minimal(base_size = 10) +
    labs(x = "mean |(SHAP value)|") +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 0),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.line.x = element_line(color = 'black'),
          axis.ticks.x = element_line(color = 'black'))

  egg::ggarrange(hm, bp1, bp2,
                 nrow = 1,
                 widths = c(1, 3, 3))
}
```

``` r
generate_heatmap_and_barplots(Phylum_data)
```

<figure>
<img src="Marker_Taxa_files/figure-gfm/Phylum%20hm-1.png"
alt="Phylum Enrichment, Significance and Importance by DAA Tools and SHAP Values of the Grass-Drought Dataset. Binary heatmap showing the enrichment of the top 15 significant phyla from ANCOM-BC2 between Control (blue) and Drought (red) groups for the three methods used for DAA (DESeq2, ANCOM-BC2, ALDEx2) with an alpha &lt;0.05, and SHAP values obtained from the RFC. Corresponding bar plots comparing -log_{10}(p_adjust) values (orange) and mean(|SHAP values|) (green)." />
<figcaption aria-hidden="true"><strong>Phylum Enrichment, Significance
and Importance by DAA Tools and SHAP Values.</strong> Binary heatmap
showing the enrichment of the top 15 significant phyla from ANCOM-BC2
between Control (blue) and Drought (red) groups for the three methods
used for DAA (DESeq2, ANCOM-BC2, ALDEx2) with an alpha &lt;0.05, and
SHAP values obtained from the RFC. Corresponding bar plots comparing
<span
class="math inline"> − <em>l</em><em>o</em><em>g</em><sub>10</sub></span>(p_adjust)
values (orange) and mean(|SHAP values|) (green).</figcaption>
</figure>

``` r
generate_heatmap_and_barplots(Class_data)
```

<figure>
<img src="Marker_Taxa_files/figure-gfm/Class%20hm-1.png"
alt="Class Enrichment, Significance and Importance by DAA Tools and SHAP Values of the Grass-Drought Dataset. Binary heatmap showing the enrichment of the top 25 significant classes from ANCOM-BC2 between Control (blue) and Drought (red) groups for the three methods used for DAA (DESeq2, ANCOM-BC2, ALDEx2) with an alpha &lt;0.05, and SHAP values obtained from the RFC. Corresponding bar plots comparing -log_{10}(p_adjust) values (orange) and mean(|SHAP values|) (green)." />
<figcaption aria-hidden="true"><strong>Class Enrichment, Significance
and Importance by DAA Tools and SHAP Values.</strong> Binary heatmap
showing the enrichment of the top 25 significant classes from ANCOM-BC2
between Control (blue) and Drought (red) groups for the three methods
used for DAA (DESeq2, ANCOM-BC2, ALDEx2) with an alpha &lt;0.05, and
SHAP values obtained from the RFC. Corresponding bar plots comparing
<span
class="math inline"> − <em>l</em><em>o</em><em>g</em><sub>10</sub></span>(p_adjust)
values (orange) and mean(|SHAP values|) (green).</figcaption>
</figure>

``` r
generate_heatmap_and_barplots(Order_data)
```

<figure>
<img src="Marker_Taxa_files/figure-gfm/Order%20hm-1.png"
alt="Order Enrichment, Significance and Importance by DAA Tools and SHAP Values of the Grass-Drought Dataset. Binary heatmap showing the enrichment of the top 25 significant orders from ANCOM-BC2 between Control (blue) and Drought (red) groups for the three methods used for DAA (DESeq2, ANCOM-BC2, ALDEx2) with an alpha &lt;0.05, and SHAP values obtained from the RFC. Corresponding bar plots comparing -log_{10}(p_adjust) values (orange) and mean(|SHAP values|) (green)." />
<figcaption aria-hidden="true"><strong>Order Enrichment, Significance
and Importance by DAA Tools and SHAP Values.</strong> Binary heatmap
showing the enrichment of the top 25 significant orders from ANCOM-BC2
between Control (blue) and Drought (red) groups for the three methods
used for DAA (DESeq2, ANCOM-BC2, ALDEx2) with an alpha &lt;0.05, and
SHAP values obtained from the RFC. Corresponding bar plots comparing
<span
class="math inline"> − <em>l</em><em>o</em><em>g</em><sub>10</sub></span>(p_adjust)
values (orange) and mean(|SHAP values|) (green).</figcaption>
</figure>

``` r
generate_heatmap_and_barplots(Family_data)
```

<figure>
<img src="Marker_Taxa_files/figure-gfm/Family%20hm-1.png"
alt="Family Enrichment, Significance and Importance by DAA Tools and SHAP Values of the Grass-Drought Dataset. Binary heatmap showing the enrichment of the top 25 significant families from ANCOM-BC2 between Control (blue) and Drought (red) groups for the three methods used for DAA (DESeq2, ANCOM-BC2, ALDEx2) with an alpha &lt;0.05, and SHAP values obtained from the RFC. Corresponding bar plots comparing -log_{10}(p_adjust) values (orange) and mean(|SHAP values|) (green)." />
<figcaption aria-hidden="true"><strong>Family Enrichment, Significance
and Importance by DAA Tools and SHAP Values.</strong> Binary heatmap
showing the enrichment of the top 25 significant families from ANCOM-BC2
between Control (blue) and Drought (red) groups for the three methods
used for DAA (DESeq2, ANCOM-BC2, ALDEx2) with an alpha &lt;0.05, and
SHAP values obtained from the RFC. Corresponding bar plots comparing
<span
class="math inline"> − <em>l</em><em>o</em><em>g</em><sub>10</sub></span>(p_adjust)
values (orange) and mean(|SHAP values|) (green).</figcaption>
</figure>

``` r
generate_heatmap_and_barplots(Genus_data)
```

<figure>
<img src="Marker_Taxa_files/figure-gfm/Genus%20hm-1.png"
alt="Genera Enrichment, Significance and Importance by DAA Tools and SHAP Values of the Grass-Drought Dataset. Binary heatmap showing the enrichment of the top 25 significant genera from ANCOM-BC2 between Control (blue) and Drought (red) groups for the three methods used for DAA (DESeq2, ANCOM-BC2, ALDEx2) with an alpha &lt;0.05, and SHAP values obtained from the RFC. Corresponding bar plots comparing -log_{10}(p_adjust) values (orange) and mean(|SHAP values|) (green)." />
<figcaption aria-hidden="true"><strong>Genera Enrichment, Significance
and Importance by DAA Tools and SHAP Values.</strong> Binary heatmap
showing the enrichment of the top 25 significant genera from ANCOM-BC2
between Control (blue) and Drought (red) groups for the three methods
used for DAA (DESeq2, ANCOM-BC2, ALDEx2) with an alpha &lt;0.05, and
SHAP values obtained from the RFC. Corresponding bar plots comparing
<span
class="math inline"> − <em>l</em><em>o</em><em>g</em><sub>10</sub></span>(p_adjust)
values (orange) and mean(|SHAP values|) (green).</figcaption>
</figure>
