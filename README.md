### Metabolomics Analysis Pipeline

#### **Overview**
This repository contains an end-to-end computational pipeline for metabolomics data analysis, integrating **statistical processing**, **biomarker discovery**, and **survival modeling**. Designed for studies of disease progression (e.g., cancer), it processes raw metabolomic profiles into clinically interpretable insights.

---

### **Pipeline Workflow**
#### 1. **MetaboAnalyst Processing** (`metaboanalyst.R`)
_Input_: Raw metabolomics data (`data.txt`)  
_Output_: Normalized data, statistical models, and visualizations  

**Key Steps**:
```mermaid
graph LR
A[Raw Data] --> B[Preprocessing]
B --> C[Normalization<br>MedianNorm + LogNorm]
C --> D[Univariate Analysis<br>Fold Change + T-tests]
D --> E[PCA/PLS-DA<br>Dimensionality Reduction]
E --> F[ROC Analysis<br>Biomarker Validation]
F --> G[Results Export]
```

**Output Structure**:
```
results/
├── 1_normalisation/      # QC plots & normalized data
├── 2_univariate/         # Fold change & t-test results
├── 3_PCA/                # PCA scree/loading plots
├── 4_PLS/                # PLS model performance
├── 5_ROC/                # ROC curves & biomarker importance
└── data_filtered.csv     # Significant metabolites
```

---

#### 2. **Survival Analysis** (`survival.R`)
_Input_: Normalized metabolomics data + clinical metadata  
_Output_: Survival models linking metabolites to patient outcomes  

**Key Analyses**:
- **Univariate Screening**:
  - Cox PH models for continuous variables (age, AFP)
  - Log-rank tests for categorical variables (etiology, BCLC stage)
- **Survival Modeling**:
  - Overall Survival (OS) and Progression-Free Survival (PFS) analysis
  - Hazard ratios (HR) for high/low metabolite expression
  - FDR-adjusted p-values (`p.adj < 0.25`)

**Outputs**:
- `univariate.xlsx`: Clinical covariate associations  
- `OS_no_adjustments.xlsx`/`PFS_no_adjustments.xlsx`: Metabolite survival signatures  
- Kaplan-Meier plots for significant metabolites  

---

### **Dependencies**
| **Package**      | Purpose                          |
|------------------|----------------------------------|
| `MetaboAnalystR` | Core metabolomics processing     |
| `survival`       | Survival modeling                |
| `survminer`      | Survival visualization           |
| `tidyverse`      | Data manipulation                |
| `openxlsx`       | Excel report generation          |

---

### **Troubleshooting**
- **Rserve Conflicts**: Scripts include `system("kill -9 PID")` to resolve port conflicts  

---

### **References**
1. Pang Z, Xu L, Viau C, Lu Y, Salavati R, Basu N, Xia J. MetaboAnalystR 4.0: a unified LC-MS workflow for global metabolomics. Nat Commun. 2024 May 1;15(1):3675. doi: 10.1038/s41467-024-48009-6. PMID: 38693118; PMCID: PMC11063062.





