rm(list=ls())
require(marray)
require(survival)
require(stringr)
require(gplots)
require(marray) 
require(ConsensusClusterPlus)
require(plyr)
require(ggplot2)
require(multtest)
require(estimate)
require(psych)
require(writexl)
if(!require(e1071)) install.packages('e1071', dependencies = TRUE)
if(!require(caret)) install.packages('caret', repos = "http://cran.us.r-project.org")
if(!require(dplyr)) install.packages('dplyr', repos = "http://cran.us.r-project.org")
if(!require(tidyverse)) install.packages('tidyverse', repos = "http://cran.us.r-project.org")
if(!require(readxl)) install.packages('readxl', repos = "http://cran.us.r-project.org")
if(!require(xlsx)) install.packages('xlsx', repos = "http://cran.us.r-project.org")
if(!require(glmnet)) install.packages('glmnet', repos = "http://cran.us.r-project.org")
if(!require(dummies)) install.packages('dummies', repos = "http://cran.us.r-project.org")
if(!require(mRMRe)) install.packages('mRMRe', repos = "http://cran.us.r-project.org")
if(!require(geometry)) install.packages('geometry', repos = "http://cran.us.r-project.org")
if(!require(survival)) install.packages('survival', repos = "http://cran.us.r-project.org")
if(!require(survminer)) install.packages('survminer', repos = "http://cran.us.r-project.org")
if(!require(ggfortify)) install.packages('ggfortify', repos = "http://cran.us.r-project.org")
if(!require(ggpubr)) install.packages('ggpubr', repos = "http://cran.us.r-project.org")
if(!require(MatchIt)) install.packages('MatchIt', repos = "http://cran.us.r-project.org")
if(!require(writexl)) install.packages('writexl', repos = "http://cran.us.r-project.org")
if(!require(ggcorrplot)) install.packages('ggcorrplot', repos = "http://cran.us.r-project.org")


local_path_GSEA = '/Users/rinading/Desktop/UCLA/WQE/data/GSEA/'
genomic = read_excel(paste0(local_path_GSEA, "matched_luad_genomic_all.xlsx"), col_names = TRUE)
risk_scores = read_excel(paste0(local_path_GSEA, "luad_risk_scores_all.xlsx"), col_names = TRUE)
geneids = read_excel(paste0(local_path_GSEA, "genesymbols.xlsx"), col_names = TRUE)

#=============================
# Pre select a set of genes #
#=============================

spearman_result = data.frame(matrix(nrow = nrow(genomic), ncol = 2))
colnames(spearman_result) = c("corr", "p")
for (gene_index in 1:nrow(genomic)){
  result = corr.test(risk_scores, t(genomic[gene_index, ]), method = "spearman")
  spearman_result[gene_index, 1] = result$r
  spearman_result[gene_index, 2] = result$p
}

genes_to_keep = c()
for (result_index in 1:nrow(spearman_result)){
  if (spearman_result[result_index, 2] < 0.05){
    genes_to_keep = c(genes_to_keep, result_index)
  }
}

geneids_selected = geneids[genes_to_keep, ]
write_xlsx(geneids_selected, paste0(local_path_GSEA, "luad_selected_genes.xlsx"))

#=================================
# Run the GSEA library to get the ES #
#=================================
setwd('/Users/rinading/Desktop/CCIPD/immune_project/GEA/ssGSEAProjection/')
source('/Users/rinading/Desktop/CCIPD/immune_project/GEA/ssGSEAProjection/ssGSEAProjection.R')
source('/Users/rinading/Desktop/CCIPD/immune_project/GEA/ssGSEAProjection/ssGSEAProjection.Library.R')
source('/Users/rinading/Desktop/CCIPD/immune_project/GEA/ssGSEAProjection/common.R')

ssGSEA.project.dataset(javaexec = "ssgseaprojection.jar", jardir = "/Users/rinading/Desktop/CCIPD/immune_project/GEA/",
                       input.ds = "/Users/rinading/Desktop/UCLA/WQE/data/GSEA/INPUT_LUAD.gct",
                       output.ds = "/Users/rinading/Desktop/UCLA/WQE/data/GSEA/OUTPUT_LUAD.gct", gene.sets.dbfile.list = "/Users/rinading/Desktop/UCLA/WQE/data/GSEA/SIGNATURE_LUAD.gmx")

#=================================
# Correlate ES with radiomic #
#=================================

enrichment_scores = read_excel("/Users/rinading/Desktop/UCLA/WQE/data/enrichment_radiomic_corr/enrichment_score_luad.xlsx", col_names = TRUE)
#enrichment_scores_values = enrichment_scores[, 2:94]
enrichment_scores_values = t(enrichment_scores_values)
radiomic = read_excel("/Users/rinading/Desktop/UCLA/WQE/data/enrichment_radiomic_corr/selected_radiomic_luad.xlsx", col_names = TRUE)
radiogenomic_result_corr = data.frame(matrix(nrow = ncol(radiomic), ncol = ncol(enrichment_scores_values)))
radiogenomic_result_p = data.frame(matrix(nrow = ncol(radiomic), ncol = ncol(enrichment_scores_values)))

for (radiomic_index in 1:ncol(radiomic)){
  for (genomic_index in 1:ncol(enrichment_scores_values)){
    result = corr.test(radiomic[, radiomic_index], enrichment_scores_values[, genomic_index], method = "spearman")
    radiogenomic_result_corr[radiomic_index, genomic_index] = result$r
    radiogenomic_result_p[radiomic_index, genomic_index] = result$p
    
  }
}

radiogenomic_result_p_corrected =
  radiogenomic_result_p %>%
  as.matrix %>%
  p.adjust(method = 'fdr') %>%
  matrix(nrow = ncol(radiomic), ncol = nrow(enrichment_scores))

rownames(radiogenomic_result_corr) = colnames(radiomic)
colnames(radiogenomic_result_corr) = unlist(enrichment_scores[, 1])
rownames(radiogenomic_result_p_corrected) = colnames(radiomic)
colnames(radiogenomic_result_p_corrected) = unlist(enrichment_scores[, 1])
rownames(radiogenomic_result_p) = colnames(radiomic)
colnames(radiogenomic_result_p) = unlist(enrichment_scores[, 1])
# Margin Order: top, right, bottom, left
ggcorrplot(as.matrix(radiogenomic_result_corr), p.mat = as.matrix(radiogenomic_result_p_corrected), insig = "blank", lab = TRUE, lab_size = 2.5) + theme(axis.text.y=element_text(size=9), axis.text.x=element_text(size=9))
ggsave(filename = "/Users/rinading/Desktop/luad.png", device='tiff', dpi=700, bg = "#FFFFFF")

