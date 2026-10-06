#' Power Calculation from gene expression data information
#'
#' This function takes the required input information such as count data, sample data, etc. to calculate the power.
#' It filters the input count data, performs DESeq2 analysis to calculate differentially expressed genes (DEGs), and
#' then calculates the power of detecting these DEGs based on simulations.
#'
#' @param countDat A matrix or data frame of raw count data where rows represent genes and columns represent samples.
#' @param smplDat A data frame of sample information, with at least a `condition` column that specifies the experimental condition of each sample.
#' @param alpha The significance level (FDR threshold) used to identify differentially expressed genes. Default is 0.05.
#' @param thrsholdFC The threshold for the absolute value of log2 fold change used to filter DEGs. Default is 2.
#' @param inptNoOfReplicates The input number of replicates based on which the power will be calculated. Default is 3.
#' @param sims The number of simulations to run for power calculation. Default is 10.
#'
#' @return A data frame containing the calculated power values and related parameters.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results counts dispersions
#' @importFrom ssizeRNA check.power
#' @importFrom dplyr filter
#' @importFrom stats complete.cases
#' @importFrom utils write.table
#' @importFrom stats na.omit
#' @importFrom Rdpack reprompt
#'
#' @references
#' Bi et al. (2016) \doi{10.1186/s12859-016-0994-9}
#' Love et al. (2014) \doi{10.1186/s13059-014-0550-8}
#'
#'
#' @details
#' Example files included with this package:
#' - `exmplCountDat.csv`: A toy dataset with count data.
#' - `exmplSampleDat.csv`: A sample dataset with metadata.
#'
#' These files are stored in the `inst/extdata` directory and can be accessed
#' using the `system.file()` function in R.
#'
#' @examples
#' \donttest{
#' # Load example files
#' countDatPath <- system.file("extdata", "exmplCountDat.csv", package = "HEssRNA")
#' smplDatPath <- system.file("extdata", "exmplSampleDat.csv", package = "HEssRNA")
#'
#' if (file.exists(countDatPath) && file.exists(smplDatPath)) {
#'   countDat <- read.csv(countDatPath)
#'   smplDat <- read.csv(smplDatPath)
#'
#'   result <- powerCalc(countDat, smplDat)
#'   print(result$PowerResults)
#' } else {
#'   warning("Example data files not found.")
#' }
#' }
#'
#' @export
powerCalc <- function(countDat, smplDat, alpha = 0.05, thrsholdFC = 2, inptNoOfReplicates = 3, sims = 10)
{
  rownames(smplDat) <- smplDat[,1]
  smplDat <- smplDat[,-1]
  smplDat$condition <- as.factor(smplDat$condition)
  rownames(countDat) <- countDat[,1]
  countDat <- countDat[,-1]
  #filter genes with 0 counts in all samples
  countDatFltrd <- countDat[rowSums(countDat) > 0, ]
  NumOfFilteredGenes <- nrow(countDatFltrd)
  #Run DESeq2 for DEGs calculation
  dds <- DESeqDataSetFromMatrix(countData = countDatFltrd , colData = smplDat, design = ~condition)
  dds <- DESeq(dds)
  # Extract results and convert to data frame for consistency
  res <- as.data.frame(na.omit(results(dds)))
  # The idea is to remove NA value from pvalue column which may cause in unnecessary error such as smooth.spline
  res <- res[complete.cases(res$pvalue),]
  res <- res[complete.cases(res$padj),]
  # Filter DEGs based on FDR (alpha) and log2FC threshold
  DEGsFltrd <- res[res$padj < alpha & abs(res$log2FoldChange) >= thrsholdFC, ]
  DEGsFltrd <- as.data.frame(DEGsFltrd)  # Ensure it's a data frame
  noOfDEGs <- nrow(DEGsFltrd) #Nuber of DEGs after applying both filters
  message("Number of DEGs:", noOfDEGs, "\n")
  # Identify upregulated genes (log2FoldChange >= thrsholdFC)
  upRegDEGs <- DEGsFltrd[DEGsFltrd$log2FoldChange >= thrsholdFC, ]
  noOfUpRegDEGs <- nrow(upRegDEGs)
  # Identify downregulated genes (log2FoldChange <= -thrsholdFC)
  downRegDEGs <- DEGsFltrd[DEGsFltrd$log2FoldChange <= -thrsholdFC, ]
  noOfDownRegDEGs <- nrow(downRegDEGs)
  # Output the number of upregulated and downregulated genes
  message("Number of Upregulated DEGs:", noOfUpRegDEGs, "\n")
  message("Number of Downregulated DEGs:", noOfDownRegDEGs, "\n")
  #inptPwr calculation
  #calculate gene-wise mean and dispersion using DESeq2
  #Obtain average read count (normalized):
  norm_counts <- counts(dds, normalized=TRUE)
  #Obtain gene-wise dispersions:
  dispersions1 <- dispersions(dds)
  #calculate mean for normalized counts
  mu <- rowMeans(norm_counts)
  nGenes <- nrow(countDatFltrd)
  pi0 = (nGenes-noOfDEGs)/nGenes
  up = (noOfDEGs - noOfUpRegDEGs)/noOfDEGs
  if(is.nan(up)){
    message("There is no significant up-regulated gene")
  }else{
    message("pi0: ", pi0)
    message("nGenes:", nGenes, "\n")
    message("noOfDEGs:", noOfDEGs)
    message("Calculating the inptPwr...\n")
    pwrCalculated.list <- check.power(nGenes = nGenes, pi0 = pi0, m = inptNoOfReplicates, mu = mu, disp = dispersions1 , fc = thrsholdFC, up = up,
                                      replace = TRUE, fdr = alpha, sims)
    # Convert list to dataframe without row names
    pwrCalculated.df <- data.frame(
      Parameter = names(pwrCalculated.list),
      Value = unlist(pwrCalculated.list),  # Use unlist to flatten the list to a vector
      stringsAsFactors = FALSE,  # Avoid converting strings to factors
      row.names = NULL  # Set row names to NULL
    )
    # Combine results into a list and return
    res.lst <- list(
      PowerResults = pwrCalculated.df,
      DESeqResults = res,
      FilteredDEGs = DEGsFltrd
    )
    return(res.lst)
  }
}
