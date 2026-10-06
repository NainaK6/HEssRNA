#' Process Data Frame in In-House Format for Model Building
#'
#' This function takes a data frame in an in-house format and processes it
#' to make it in longer format and round the value of the power to 3 digits
#' for building a model. It reshapes the data from a wide format to a long
#' format, extracting and manipulating columns related to replicate numbers and power values.
#' This function is needed when user has a data frame similar to the in-house format.
#' For the purpose of creating model the user should also have Heritability class and log fold change value too.
#'
#' @param df4modelInhouseFmt A data frame containing the input data in in-house format.
#'        The columns should include replicate columns named starting with "R" (e.g., R1, R2, etc.).
#'
#' @return A data frame in long format with columns:
#'   \item{NoOfReplicates}{Numeric representation of the replicate number extracted from column names (R1, R2, etc.).}
#'   \item{pwr}{Power values rounded to 3 decimal places corresponding to the replicate number.}
#'
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr starts_with
#'
#'
#' @examples
#' # Example of usage:
#' df <- data.frame(
#'   Gene = c("Gene1", "Gene2"),
#'   R1 = c(0.85, 0.90),
#'   R2 = c(0.88, 0.91),
#'   R3 = c(0.83, 0.89)
#' )
#' result <- prcesDF4modelInhouse(df)
#' print(result)
#'
#' @export
prcesDF4modelInhouse <- function(df4modelInhouseFmt)
{
  dfLong <- df4modelInhouseFmt %>%
    pivot_longer(cols = starts_with("R"), names_to = "NoOfReplicates", values_to = "pwr") %>%
    mutate(NoOfReplicates = as.numeric(sub("R", "", NoOfReplicates)))
  dfLong$pwr <- (round(dfLong$pwr, 3) + 0.001)
  dfLong <- as.data.frame(dfLong)
  return(dfLong)
}

