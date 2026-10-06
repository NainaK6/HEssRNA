#' Generate a Linear Model for Sample Size Prediction
#'
#' This function generates a linear regression model to predict the number of replicates (`NoOfReplicates`) based on heritability, power, fold change, and tissue type. The model is generated depending on whether the tissue information is provided in the data. The function returns the fitted model.
#'
#' @param df4model A data frame containing the input data for the model. It should include the following columns: `NoOfReplicates`, `HeritabilityValue`, `pwr`, `FoldChange`, and optionally, `Tissue`.
#'
#' @param heritabilityClass A character value indicating the class of heritability used for filtering the data.
#'
#' @param inptPwr A numeric value representing the power used in the model.
#'
#' @param fc A numeric value representing the fold change used in the model.
#'
#' @param trait An optional parameter to specify the trait. If provided, it can be used for further filtering, but it's not currently used in the function.
#'
#' @param tissue An optional parameter specifying the tissue type. If provided, the model will include the tissue information in the regression. If not provided, the model will exclude tissue information.
#'
#' @return A linear model object (`lm` class), which contains the fitted linear regression model for the number of replicates prediction.
#'
#' @importFrom stats lm
#' @importFrom Rdpack reprompt
#'
#' @references
#' Sun et al. (2017) \doi{10.1093/nar/gkx204}
#'
#' @examples
#' \donttest{
#' # Example usage:
#' df4modelInpt <- data.frame(
#'     NoOfReplicates = c(3, 5, 7, 9, 11),
#'     HeritabilityClass = c("high", "mid", "low", "high", "mid"),
#'     HeritabilityValue = c(0.5, 0.6, 0.7, 0.5, 0.6),
#'     pwr = c(0.8, 0.9, 0.85, 0.88, 0.86),
#'     FoldChange = c(2, 3, 2.5, 2.8, 3.2),
#'     Tissue = c("Liver", "Liver", "Kidney", "Liver", "Kidney")
#' )
#'
#' # Fit the model
#' model <- smplSizPredModel(
#'     df4model = df4modelInpt,
#'     heritabilityClass = "high",
#'     inptPwr = 0.8,
#'     fc = 2,
#'     tissue = "Liver"
#' )
#'
#' # Summarize the results
#' summary(model)
#' }
#'
#' @export
smplSizPredModel <- function(df4model = df4modelInpt, heritabilityClass, inptPwr, fc, trait=NULL, tissue = NULL)
{
  suppressWarnings({
    if (!is.null(tissue) & nrow(df4model[df4model$HeritabilityClass == heritabilityClass & df4model$Tissue == tissue, ]) > 0) {
      model <- lm(NoOfReplicates ~ HeritabilityValue + pwr + FoldChange + Tissue, data = df4model)
    } else {
      model <- lm(NoOfReplicates ~ HeritabilityValue + pwr + FoldChange, data = df4model)
    }
  })
  return(model)
}
