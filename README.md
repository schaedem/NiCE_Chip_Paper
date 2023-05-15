# NiCE_Chip_Paper
Contains data and code to generate figures/analysis for Schaedel et al. (2023) in _PLOS ONE_: <i> Temporal assessment of N-cycle microbial functions in a tropical agricultural soil using gene co-occurrence networks </i>

## Data Files & Data Wrangling

### 16S_data.csv: 
16S rRNA abundance as assessed by conventional qPCR

### final_nice_std_data.csv
processed data from HT-qPCR. Raw gene abundances (conc) were transformed to log10 (log_conc) and used to obtain abundances normalized to g soil (log_conc_truesoil) for subsequent analyses

### assay_list_2.csv: 
metadata for HT-qPCR assays, including target gene/organism (gene_org), enzymatic pathway (pathway), color for graphing (colorcode), and abbreviation for manusript (acronym)

### nice_matrix.Rdata:
Matrix format of normalized gene abundances (log_conc_truesoil), pre-requisite for Spearman dissimilarity matrix (R data object)

### dist_spearman.Rdata:
Spearman dissimilarity matrix (all timepoints, both locations) used to construct networks and perform multivariate analyses (R data object)

### correlation_dataframes.R:
Code to make dataframes for each location x timepoint which included all pairwise Spearman correlations, including rho, P, and sample metadata 

### full_dat.Rdata:
final data used in analysis plus metadata (R data object)

### networks_and_statistics.R:
compile network statistics into dfs for easy plotting

## Analysis Files

### permanova_community_differences.R:
Significance testing for differences in N cycle gene structure - pre-requisite for co-occurrence analysis

### intersection_networks.R:
Create graphs for each loc x timepoint, then find the interesection of those graphs. Calculate network statistics for each intersection graph. 

### erdos_renyi_random_graphs.R:
random network statistics for Table 2

### permanova_community_differences.R

### cooccurrence_permanova_method2.R:
Analysis for Fig 2

### supplemental_network_stats.R:
network statistics for supplemental table 4

## Figures

### cooccurrence_permanova_method2.R:
Fig 2

### coccurrence_dbRDA_2.R:
Fig 3

### network_graphs.R: 
Base functions to generate network graphs for Fig 4

### intersection_networks.R:
Generate individual network graphs and collate for Fig 4

### intersection_graph_betweenness_closeness.R:
Fig 5

### correlogram_heatmap.R:
Supplemental Fig 5

### timepoint_cocurrence_dbRDA.R:
Supplemental Fig 6


