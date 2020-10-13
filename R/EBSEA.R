#' @title Exon Based Startegy for Expression Analysis of genes
#' @description EBSEA takes the filtered raw exon-level read counts
#' as input, normalizes and performs a two-group statistical
#' comparison with DESeq2. The exon-level results are
#' aggregated to the gene-level using empirical Brown’s method.
#' The samples in the two groups can be paired.
#' @param data A dataframe of raw exon-counts
#' @param columnData A dataframe indicated the groups of the samples.
#' @param design Design matrix (see more information od design
#' matrixes in DESeq2 reference manual)
#' @param test The statistical test to be carried out. It can be
#' either Wald or Likelihood Ratio Test. For further details about
#' the methods you can look into DESeq2 refernce manual. Default: Wald
#' @param contrasts a character vector with exactly three elements:
#' the name of a factor in the design formula, the name of the
#' numerator level for the fold change, and the name of the
#' denominator level for the fold change Default: NULL
#' @param plot A logical value indicating a volcano plot
#' is produced. Default: FALSE
#' @return The function returns a list containing containing
#' exon and gene-level results. ExonTable is a data frame
#' containing an average expression, log2 fold-change,
#' p-value and adjusted p-value. GeneTable is a data frame
#' containing the gene-level p-value, and adjusted-value.
#' Other returned elements include the raw and normalised
#' exon-level read counts, group information and design matrix used.
#' @references Laiho, A., & Elo, L. L. (2014). A note on an
#' exon-based strategy to identify differentially expressed
#' genes in RNA-seq experiments. PloS One, 9(12), e115964.
#' @examples # The exon-based analysis for unpaired samples can be performed as follows:
#' data(exonCounts)
#' group <- data.frame('group' = as.factor(c('G1', 'G1', 'G1', 'G2', 'G2', 'G2', 'G2')))
#' row.names(group) <- colnames(exonCounts)
#' design <- ~group
#' ebsea.out <- EBSEA(exonCounts, group, design)
#' @examples # The exon-based analysis for paired samples with contrast provided can be performed as follows:
#' data(exonCounts)
#' group <- data.frame('group' = as.factor(c('G1', 'G1', 'G1', 'G2', 'G2', 'G2', 'G2')),
#'  'paired' = as.factor(c(1,2,3,1,2,3,3)))
#' row.names(group) <- colnames(exonCounts)
#' design <- ~group
#' contrastInfo <- c('group', 'G2', 'G1')
#' ebsea.out <- EBSEA(exonCounts, group, design, contrasts = contrastInfo)
#' @import DESeq2
#' @import EmpiricalBrownsMethod
#' @importFrom stats p.adjust
#' @export EBSEA

EBSEA <- function(data, columnData, design, test = 'Wald', contrasts = NULL, plot = FALSE){
  message("Checking Parameters")
  if(is.null(data)) {
    stop('No count data found..')
  }
  # Checking the the design
  
  message("Performing Statistical testing of Exons")
  DESeqData <- DESeqDataSetFromMatrix(data, columnData, design)
  DESeqData <- estimateSizeFactors(DESeqData)
  norm.counts <- counts(DESeqData, normalized = TRUE)
  DESeqData <- estimateDispersions(DESeqData)
  if(test == 'Wald'){
    message('Using Wald test for testing')
    DESeqData <- nbinomWaldTest(DESeqData, betaPrior = FALSE)
  } else if(test == 'LRT'){
    message('Using LRT test for testing')
    DESeqData <- nbinomLRT(DESeqData, reduced = ~ 1)
  } else{
    message('No test provided')
  }
  if(is.null(contrasts)){
    exon.table <- as.data.frame(results(DESeqData))
  }else{
    exon.table <- as.data.frame(results(DESeqData, contrast = contrasts))
  }
  
  # Aggreagting to gene level results
  message('Aggregating to Gene-level results')
  exon.table$gene_id <- gsub(':.*', '', row.names(exon.table))
  exon.table <- exon.table[order(row.names(exon.table)), ]
  norm.counts <- norm.counts[order(row.names(norm.counts)), ]
  norm.counts.stats <- data.frame(norm.counts,
                                  'P.Value' = exon.table$pvalue,
                                  'gene_id' = exon.table$gene_id)
  norm.counts.stats <- norm.counts.stats[!is.na(norm.counts.stats$P.Value), ]
  gene.table <- do.call(rbind, by(norm.counts.stats,
                                  list(norm.counts.stats$gene_id), calculate_ebm))
  gene.table <- data.frame('Gene' = row.names(gene.table),
                           'pvalue' = gene.table[,'P_test', drop= FALSE],
                           'padj' = p.adjust(gene.table[,'P_test', drop= FALSE]))
  
  ## editing the genetable
  
  
  result <- list("ExonTable" = exon.table,
                 "GeneTable" = gene.table,
                 "RawCounts" = data,
                 "NormCounts" = norm.counts,
                 "Group" = columnData$group,
                 "design" = design)
  
  ## MA plot for exon-level results
  
  ## Volcano plot for gene-level results
  
  tryCatch(if(plot == TRUE) {
    with(result$ExonTable, plot(result$ExonTable$log2FoldChange, -log2(result$ExonTable$pvalue),
                                pch = 20, main = "Volcano plot"))
    abline(v = c(1, -1), lty = 3, lwd = 1)
  }, error = function(e){print(e)})
  class(result) <- "EBSEA"
  rm(exon.table)
  rm(gene.table)
  rm(norm.counts)
  rm(data)
  gc()
  message('Done')
  return(result)
}

calculate_ebm <- function(x){
  empiricalBrownsMethod(data_matrix = as.matrix(x[,seq(1, (ncol(x)-2))]),
                        p_values=x$P.Value, extra_info = TRUE)
}

