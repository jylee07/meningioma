 # Meningioma snRNA-seq Analysis Code
\
This repository contains code used for the analysis in the study:

Single-cell analysis reveals a longitudinal trajectory of meningioma evolution and heterogeneity

\
## 1. Pseudo-bulk DEG Analysis (Fig 2)\
01_DEG_analysis/fig2_pseudobulk_deseq2.R\
- Performs DESeq2-based DEG analysis and pathway enrichment.\
\
## 2. MG4 and Risk Score (Fig 2)\
\
02_Scoring/fig2_mg4_score.R \'97 ssGSEA for MG4 meta-modules.\
02_Scoring/fig3_risk_score.R \'97 Risk score from gene sets.\
\
## 3. CellChat Signaling (Fig 4)\
03_CellChat/fig4_cellchat_primary_vs_recurrent.R\
- Infers intercellular communication & compares across samples.\
\
## 4. Differential Abundance with miloR (Fig 5)\
04_miloR/fig5_milor_analysis.R\
- Tests local differential abundance using `miloR`.\
\

## For questions or support, contact:\
ghsl0102@korea.ac.kr
