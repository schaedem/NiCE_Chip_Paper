# NiCE_Chip_Paper
Contains data and code to generate figures/analysis for Schaedel et al. (in review): <i> Temporal assessment of N-cycle microbial functions in a tropical agricultural soil using gene co-occurrence networks </i>

## Data Files & Data Wrangling

### 16S_data.csv: 
16S rRNA abundance as assessed by conventional qPCR


### final_nice_std_data.csv
processed data from HT-qPCR. Raw gene abundances (conc) were transformed to log10 (log_conc) and used to obtain abundances normalized to g soil (log_conc_truesoil) for subsequent analyses

### assay_list_2.csv: 
metadata for HT-qPCR assays, including target gene/organism (gene_org), enzymatic pathway (pathway), color for graphing (colorcode), and abbreviation for manusript (acronym)

### dist_spearman.Rdata:
Spearman dissimilarity matrix (all timepoints, both locations) used to construct networks and perform multivariate analyses

### correlation_dataframes.R:
Code to make dataframes for each location x timepoint which included all pairwise Spearman correlations, including rho, P, and sample metadata 

## Analysis Files

### erdos_renyi_random_graphs.R:
random network statistics for Table 2

### permanova_community_differences.R

### cooccurrence_permanova_method2.R

## Figures

### coccurrence_dbRDA_2.R:
Fig 3

### correlogram_heatmap.R:
Supplemental Fig 5
