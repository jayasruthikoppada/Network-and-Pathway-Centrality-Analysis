# Network-and-Pathway-Centrality-Analysis
Integrative Genomic Profiling of COVID-19: A Network and Pathway Centrality Analysis

# Project Overview

The COVID-19 pandemic caused by SARS-CoV-2 has presented a need to understand its molecular mechanisms and the corresponding human immune response. This project leverages high-throughput RNA-seq data from the [GEO dataset GSE152418](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152418) to uncover the key molecular drivers, pathways, and network interactions involved in the disease. By combining differential gene expression analysis, Weighted Gene Co-expression Network Analysis (WGCNA), and pathway enrichment techniques, this study highlights the complex biological processes and potential biomarkers linked to COVID-19 severity.


# Objectives

1. Identify Molecular Mechanisms:
   - Discover differentially expressed genes (DEGs) and associated pathways in COVID-19.
2. Characterize Gene Interactions:
   - Explore gene co-expression networks to identify critical modules associated with disease severity.
3. Propose Therapeutic Targets:
   - Highlight key genes and pathways for potential interventions.
4. Facilitate Reproducibility:
   - Provide a comprehensive workflow for genomic analysis using open-source tools.

# Features
# Data Analysis
Differential Expression Analysis:
  - Identified 1,062 DEGs (774 upregulated, 288 downregulated) using DESeq2.
  - Key DEGs: `ISG15`, `IFITM3` (upregulated), and `DUSP1`, `EIF2AK2` (downregulated).
Network Analysis:
  - Created gene interaction networks using Cytoscape and ReactomeFI.
  - Centrality measures (degree, betweenness, closeness) identified hub genes.
WGCNA:
  - Uncovered gene modules correlated with clinical traits (e.g., ICU severity).
  - Identified modules like MEmagenta and MEblue as contributors to severe disease.

# Biological Insights
Pathway Enrichment:
  - Enriched pathways include cell cycle regulation, apoptosis, and cytokine signaling.
  - p53 signaling pathway and oocyte meiosis highlight stress and immune responses.
Functional Enrichment:
  - GO analysis revealed significant processes like mitotic nuclear division and chromosome segregation.

# Visualizations
- PCA plots for dimensionality reduction and variance analysis.
- Volcano plots to highlight significant DEGs.
- Cytoscape network diagrams of top DEGs.
- Gene Ontology and pathway bar graphs.

# Dataset
Source: [GEO Dataset GSE152418](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152418)
Samples: Peripheral blood mononuclear cells (PBMCs) from 34 participants (17 COVID-19 patients, 17 healthy controls).
Data Highlights:
  - Prolonged activation of plasma blasts and effector T cells.
  - Alterations in myeloid cell function and interferon-stimulated genes.

Workflow and Methodology

# Tools and Libraries
1. R:
   - DESeq2 for differential expression analysis.
   - WGCNA for gene co-expression analysis.
2. Cytoscape:
   - Network visualization with the ReactomeFI plugin.
3. Python:
   - Custom scripts for data preprocessing and visualization.
4. Pathway Enrichment**:
   - SRplot and Reactome databases for GO/KEGG analyses.

# Step-by-Step Analysis
1. Data Preprocessing: Normalize and batch correct RNA-seq data.
2. Differential Expression: Identify significant genes using adjusted p-values and log2 fold changes.
3. Network Construction: Map DEGs onto functional interaction networks using Cytoscape.
                         Compute centrality measures to identify hub genes.
4. WGCNA: Cluster co-expressed genes into modules.
          Correlate modules with clinical traits like ICU and severe outcomes.
5. Pathway Enrichment: Use hypergeometric tests to link DEGs to GO terms and KEGG pathways.


# Results

# Key Findings
- Differentially Expressed Genes:
  - Significant upregulation of antiviral response genes (`ISG15`, `IFITM3`).
  - Downregulation of metabolic and cellular stress genes (`DUSP1`, `EIF2AK2`).
- **Gene Co-expression Modules**:
  - Modules like MEmagenta correlated with severe outcomes.
  - MEturquoise module potentially associated with protective responses.
- **Pathway Insights**:
  - Enrichment in cell cycle and p53 signaling pathways points to disrupted cell proliferation and immune stress responses.

 Visualizations
- Figures**:
  - PCA plots for sample variance.
  - WGCNA dendrograms and module-trait association plots.
  - Cytoscape networks highlighting central genes (`PLK1`, `KIF14`, `CDK2`).


# Repository Structure

Integrative-COVID19-Genomics/
├── data/                       # RNA-seq data files
├── analysis/                   # R scripts for DESeq2 and WGCNA
├── results/                    # Output files and processed results
├── visuals/                    # Figures, plots, and diagrams
├── docs/                       # Project documentation
├── README.md                   # Project overview


Installation:

# Prerequisites
Software:
  - R with DESeq2, WGCNA, ggplot2 installed.
  - Cytoscape with ReactomeFI plugin.
Dataset:
  - Download GSE152418 data from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152418).


# References
 [GEO Dataset GSE152418](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152418)
 [Reactome Pathway Database](https://reactome.org/)
 Key supporting literature cited in `docs/references.md`.

