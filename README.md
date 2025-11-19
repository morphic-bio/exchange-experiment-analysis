# Exchange Experiment Analysis

This repository contains data and analysis pipelines for four transcription factor (TF) knockouts studied by four KO strategies across multiple differentiation systems. The goal is to compare KO strategies and use the KO data to understand gene function.

---

## Overview of Knockout Genes and Systems

| Gene   | Center | Differentiation system                   | Time point | Likely cell types                          | Assays                          | Known target genes?             |
|--------|--------|-------------------------------------------|------------|--------------------------------------------|---------------------------------|---------------------------------|
| **EOMES** | MSK    | hPSC → DE differentiation               | Day 3      | Definitive endoderm                        | Flow cytometry, RNA-seq          | SOX17                          |
| **GCM1**  | JAX    | Extra-embryonic differentiation         | Day 6      | PrSyn cells                                | RNA-seq, fusion index            | PGF, PRR9, TWIST2, LMO2         |
| **ISL1**  | NWU    | Cardiomyocyte differentiation           | Day 32     | Cardiomyocytes, smooth muscle, fibroblast  | Impedance, Ca²⁺ imaging, RNA-seq | —                               |
| **NKX2.1**| UCSF   | Neural progenitors                     | Day 18     | ?                                          | RNA-seq                         | —                               |

---

## Common Analysis Plan

1. **Dimensionality reduction**  
   - PCA showing separation of knockouts (KOs) from controls.

2. **Comparative differential expression**  
   - Scatterplots of logFC or p-values across centers.  
   - Alternative to Venn diagrams (reservations about Venn).  



---

## Gene-Specific Analysis Plans

### EOMES
- scRNA-seq of WT and KO available; multiome of WT also available.  
- Apply gene regulatory networks (GRNs) to scRNA-seq/multiome.  
- Confirm predicted targets in bulk RNA-seq.  
- (Note: Olivier @ JAX DAV has already tested aspects of this.)

### GCM1
- scRNA-seq of WT and KO available.  
- Deconvolve bulk RNA-seq using scRNA-seq references, anticipating shifts in cell states after GCM1 KO (JAX DAV).  
- Apply GRNs to scRNA-seq and confirm targets in bulk.  
- Observation: JAX reported up-regulation of GCM1 expression after KO, but this compensation effect was not seen in other clones (e.g., UCSF, NWU).  

### ISL1
- Focus on cardiomyocyte differentiation (day 32).  
- Evaluate effects across multiple readouts (impedance, Ca^{2+} imaging, RNA-seq).  
- Cell type heterogeneity (cardiomyocytes, smooth muscle, fibroblasts) may complicate interpretation.  

### NKX2.1
- RNA-seq available for neural progenitors (day 18).  
- Cell types not fully defined — requires annotation.  

---

## Repository Structure (proposed)

```
.
├── data/              # Processed RNA-seq / scRNA-seq data
├── scripts/           # Analysis scripts (R, Python, etc.)
   ├── 01a_check_x     # Scripts to preprocess data. The x includes where the data is from (JAX, UCSF, NW, MSK) and the date the data is processed.
   ├── 02_DE_testing_x # Scripts to do DE Analysis. The x includes where the data is from (JAX, UCSF, NW, MSK) and the date the data is processed.
├── results/           # Figures, PCA, DE analysis outputs
├── docs/              # Notes, meeting summaries
└── README.md          # Project overview
```

---

## Next Steps
- Standardize QC and normalization across centers.  
- Develop shared visualization templates (PCA, scatterplots).  
- Integrate scRNA-seq/multiome GRNs with bulk results.  
- Compare compensation effects across clones and institutions.  
