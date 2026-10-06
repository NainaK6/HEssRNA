#' Predict Number of Replicates Based on Heritability, Power, and Fold Change
#'
#' This function predicts the number of replicates required for a given experiment based on heritability, power, fold change, and tissue type. The model is constructed using the provided data, and the prediction is adjusted based on the selected trait's mean heritability value. The function ensures that the predicted replicates are valid, rounding negative or unrealistic values to sensible minimums based on the heritability class.
#'
#' @param df4model A data frame containing the input data for the model. It should include the following columns: `NoOfReplicates`, `HeritabilityValue`, `pwr`, `FoldChange`, and optionally `Tissue`.
#'
#' @param hIndexMeanDFinput A data frame containing the mean heritability values for each trait. It should include at least the columns `Trait.name` and `MeanValue`.
#'
#' @param heritabilityClass A character string specifying the heritability class used for filtering and adjusting the prediction. Possible values are "low", "mid", and "high".
#'
#' @param inptPwr A numeric value representing the power used in the model.
#'
#' @param fc A numeric value representing the fold change used in the model.
#'
#' @param trait An optional parameter specifying the trait. If provided, the heritability value for the trait will be used to adjust the heritability class values.
#'
#' @param tissue An optional parameter specifying the tissue type. If provided, the model will include tissue as a factor in the regression. If not provided, tissue is excluded.
#'
#' @return A numeric value representing the predicted number of replicates. The value is rounded to the nearest whole number and adjusted to ensure it is valid for the selected heritability class.
#'
#' @importFrom stats lm na.omit predict
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
#' hIndexMeanDF <- data.frame(Trait.name = c("Trait1", "Trait2"),
#'                            MeanValue = c(0.3, 0.5))
#' NoOfReplicatesPred <- smplSizPred(df4model = df4modelInpt,
#'                       hIndexMeanDFinput = hIndexMeanDF,
#'                       heritabilityClass = "mid",
#'                       inptPwr = 0.85,
#'                       fc = 2.5,
#'                       trait = "Trait1",
#'                       tissue = "Liver")
#' print(NoOfReplicatesPred)
#' }
#'
#' @export
smplSizPred <- function(df4model = df4modelInpt, hIndexMeanDFinput = hIndexMeanDF, heritabilityClass, inptPwr, fc, trait=NULL, tissue = NULL) {

  # Step 1: Define the model and calculate heritabilityValue
  suppressWarnings({
    if (!is.null(tissue) & nrow(df4model[df4model$HeritabilityClass == heritabilityClass & df4model$Tissue == tissue, ]) > 0) {
      model <- lm(NoOfReplicates ~ HeritabilityValue + pwr + FoldChange + Tissue, data = df4model)
      heritabilityValue <- unique(df4model[(df4model$HeritabilityClass == heritabilityClass & df4model$Tissue == tissue), ]$HeritabilityValue)
    } else {
      model <- lm(NoOfReplicates ~ HeritabilityValue + pwr + FoldChange, data = df4model)
      heritabilityValue <- mean(df4model[df4model$HeritabilityClass == heritabilityClass, ]$HeritabilityValue)
    }
  })
  # Modify the value of hertiability based on the trait selected
  ##Correspond the selected trait's mean value to respective heritability class and modify the respective value. All other should be as they are.
  traitHeritabilityMeanVal <- round(hIndexMeanDFinput[hIndexMeanDFinput$Trait.name %in% trait,]$MeanValue, 3)
  traitHeritabilityMeanVal <- as.numeric(traitHeritabilityMeanVal)
  if(!is.null(trait)){
    if(traitHeritabilityMeanVal > 0.297){
      heritabilityClassValues <- c(
        "low" = 0.098,
        "mid" = 0.271,
        "high" = traitHeritabilityMeanVal
      )
    }else if(traitHeritabilityMeanVal > 0.157 & traitHeritabilityMeanVal < 0.296){
      heritabilityClassValues <- c(
        "low" = 0.098,
        "mid" = traitHeritabilityMeanVal,
        "high" = 0.447
      )
    }else if(traitHeritabilityMeanVal < 0.156){
      heritabilityClassValues <- c(
        "low" = traitHeritabilityMeanVal,
        "mid" = 0.271,
        "high" = 0.447
      )
    }
  }else{
    heritabilityClassValues <- c(
      "low" = 0.098,
      "mid" = 0.271,
      "high" = 0.447
    )
  }
  # Step 2: Apply the transformation to pwr considering both dataset HeritabilityValue and class-based values
  classBasedHeritability <- heritabilityClassValues[[heritabilityClass]]
  combinedHeritabilityEffect <- (heritabilityValue + classBasedHeritability)
  transformedPwr <- ((inptPwr/2)) * abs(log(combinedHeritabilityEffect))
  # transformedPwr <- inptPwr * abs(log(combinedHeritabilityEffect/2)) # This could also be done
  # Step 5: Create the new data frame for prediction
  if (!is.null(tissue)) {
    newDat <- data.frame(
      pwr = transformedPwr,
      HeritabilityValue = heritabilityValue,
      FoldChange = fc,
      Tissue = tissue
    )
  } else {
    newDat <- data.frame(
      pwr = transformedPwr,
      HeritabilityValue = heritabilityValue,
      FoldChange = fc
    )
  }
  # Step 3: Predict the number of replicates and suppress attributes
  suppressWarnings({
    NoOfReplicatesPred <- predict(model, newdata = newDat)
    attributes(NoOfReplicatesPred) <- NULL
  })
  # Step 4: Modify NoOfReplicatesPred so that there will be a logical value
  # Round it first
  NoOfReplicatesPred <- round(NoOfReplicatesPred)
  # In case no of replicates is in negative value modify the lowest value to be 2
  if(heritabilityClass == "high" & NoOfReplicatesPred < 2){
    NoOfReplicatesPred <- 2
  }else if(heritabilityClass == "mid" & NoOfReplicatesPred < 4){
    NoOfReplicatesPred <- 4
  }else if(heritabilityClass == "low" & NoOfReplicatesPred < 6){
    NoOfReplicatesPred <- 6
  }
  # Step 5: Return the predicted replicates rounded to the nearest whole number
  return(NoOfReplicatesPred)
}
