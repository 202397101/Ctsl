# Cathepsin L as a dual-target to mitigate muscle wasting while enhancing anti-tumor efficacy of anti-PD-L1

---

**Abstract**  
Immune checkpoint inhibitors (ICIs) have revolutionized cancer therapy; however, their use is frequently associated with immune-related adverse events (irAEs).  
In this study, anti-PD-L1 therapy exacerbates muscle wasting in tumor-bearing male mice despite its anti-tumor efficacy, accompanied by an accumulation of CD8+ T cells in muscle.  
Single-cell RNA sequencing identifies these cells as tissue-resident memory-like CD49a+ CD8+ T cells.  
While CD8+ T cell depletion prevents muscle wasting, it compromises the anti-tumor efficacy of anti-PD-L1.  
To resolve this paradox, we identify cathepsin L (CTSL) as a dual-target capable of suppressing both tumor progression and CD8+ T cell-mediated muscle wasting, through integrative transcriptomic analysis.  
Pharmacological inhibition of CTSL not only mitigates anti-PD-L1-induced muscle wasting but also further suppresses tumor growth, potentially via downregulation of BNIP3.  
Here, we show that CTSL is a dual-action target to uncouple anti-tumor efficacy from muscle-specific irAEs, offering a strategy to improve clinical outcomes of ICIs.

---

**Overview of repository**  
*Summary of R Script Execution*  
There are seven R scripts, each corresponding to a specific figure.  
To run these scripts, you must first download the necessary R data (Rdata).  
This data can be obtained from either a provided compressed file (.zip) or from the databases detailed below.  
After downloading, appropriate data processing is required before executing each script.

1. **For Figures 3B, 3E, 3F, 3G, 5E, 6A, 7A**  
To construct the protein-protein interaction (PPI) network, download the protein.links (interaction data) and protein.info (accessory data) files for *Mus musculus* from the STRING database.  
🔗 [STRING database](https://string-db.org/cgi/download?sessionId=bJCyoNzXhR2Z)

2. **For Figures 7B, 7C**  
Access the MSigDB website and download the Gene Symbols GMT file for the C2: curated gene sets.  
🔗 [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/)

3. **For Figures 3K, 3L, 6C**  
Download TCGA expression data using the TCGAbiolinks R library.  
Use the functions `getGDCprojects` and `GDCquery` with the following parameters:
```r
GDCquery(project = id, 
         data.category = "Transcriptome Profiling", 
         data.type = "Gene Expression Quantification", 
         experimental.strategy = "RNA-Seq",
         workflow.type = "STAR - Counts")

4. **For Figure 7**
From the GTEx portal, navigate to the RNA-seq section and download the following files:

RSEM_transcript_tpm.txt.gz
RSEM_transcript_expected_count.txt.gz

🔗 [GTEx]([https://www.gsea-msigdb.org/gsea/msigdb/](https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression))
