#' Calculate Mean Heritability Index for Traits
#'
#' This function processes heritability index data, filtering out empty trait names, and calculates the mean heritability for each unique trait. The resulting output is a data frame with traits and their corresponding mean heritability values.
#'
#' @param hIndexValDF A data frame containing heritability index values with at least two columns: `Trait.name` and `Heritability`.
#'                    The `Trait.name` column should contain trait identifiers, and the `Heritability` column should contain numeric heritability values.
#'
#' @return A data frame with two columns: `Trait.name` and `MeanValue`, where `MeanValue` represents the mean heritability for each trait.
#'
#' @importFrom dplyr group_by mutate ungroup
#' @importFrom magrittr %>%
#' @importFrom Rdpack reprompt
#'
#' @references
#' Hu et al. (2018) \doi{10.1093/nar/gky1084}
#'
#' @examples
#' \donttest{
#' # Example of usage:
#' hIndexValDF <- data.frame(Trait.name = c("Trait1", "Trait2", "Trait1", "Trait2"),
#'                           Heritability = c(0.5, 0.6, 0.7, 0.8))
#' result <- hIndxMeanCalc4Traits(hIndexValDF)
#' print(result)
#' }
#'
#' @export
hIndxMeanCalc4Traits <- function(hIndexValDF)
{
  hIndexVal <- hIndexValDF[!hIndexValDF$Trait.name == "",]
  hIndexVal$Heritability <- as.numeric(hIndexVal$Heritability)
  # Add a new column with mean value based on unique Group
  hIndexMean <- hIndexVal %>%
    group_by(Trait.name) %>%
    mutate(MeanValue = mean(Heritability)) %>%
    ungroup() # ungroup after the calculation
  hIndexMean <- hIndexMean[, c(1,3)]
  hIndexMean <- hIndexMean[!duplicated(hIndexMean),]
  hIndexMean <- as.data.frame(hIndexMean)
  return(hIndexMean)
}
