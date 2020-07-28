#' Predicts IDA efficacies for 2-Drug Combinations
#'
#' This function creates efficacy predictions for 2-drug combinations using monotherapy efficacy data and the assumptions of independent drug action. When data is available for multiple concentrations of each drug, efficacy predictions are made for all possible concentration combinations.
#'
#' @importFrom stats complete.cases rnorm sd aggregate
#'
#' @param Monotherapy_Data A data frame where each row contains information about the response of a single cell line to a single drug at a single concentration. Must minimally include columns containing the following information: cell line name, drug name, drug concentration, and measured drug efficacy. May optionally include a column recording the standard error (SE) of the measured drug efficacy.
#' @param Cell_Line_Name_Column A character vector of length 1 containing the name of the column in the Monotherapy_Data data frame which contains cell line names.
#' @param Drug_Name_Column A character vector of length 1 containing the name of the column in the Monotherapy_Data data frame which contains drug names.
#' @param Drug_Concentration_Column A character vector of length 1 containing the name of the column in the Monotherapy_Data data frame which contains drug concentrations.
#' @param Efficacy_Column A character vector of length 1 containing the name of the column in the Monotherapy_Data data frame which contains measured drug efficacies (i.e. percent Viability, percent Cell Growth, etc.).
#' @param LowerEfficacyIsBetterDrugEffect A logic vector of length 1 indicating whether or not lower values in Efficacy_Column indicate a more effective drug effect (i.e. for percent viability). Set TRUE if so. Otherwise, set FALSE if higher values in Efficacy_Column indicate a more effective drug response (i.e. for percent cell death).
#' @param Efficacy_Metric_Name A character vector of length 1 indicating the name of the efficacy metric being used (i.e. Percent_Viability, Percent_Growth, etc.). Used to correctly label column names in output. Defaults to "Efficacy".
#' @param Drug1 A character vector of length 1 containing the name of the first drug in the drug combination for which efficacy predictions are to be made.
#' @param Drug2 A character vector of length 1 containing the name of the second drug in the drug combination for which efficacy predictions are to be made.
#' @param Calculate_Uncertainty A logic vector of length one indicating whether or not a semi-parametric bootstrap should be performed to estimate uncertainties in the efficacy predictions based on uncertainties in the monotherapy efficacy measurements. Set TRUE if you wish to calculate uncertainties. Defaults to FALSE.
#' @param Efficacy_SE_Column A character vector of length 1 containing the name of the column in the Monotherapy_Data data frame which contains the standard errors of measured drug efficacies. Must be specified if Calculate_Uncertainty is set to TRUE.
#' @param n_Simulations A positive, integer vector of length 1 with a value >= 40 indicating the number of random samples to be drawn when calculating output efficacy prediction uncertainties. Defaults to 1000.
#' @param Calculate_IDAcomboscore_And_Hazard_Ratio A logic vector of length 1 indicating whether or not IDAcomboscores and Hazard Ratios (HRs) should be calculated between monotherapies and the drug combination. Set TRUE if so. Should only be set to TRUE for efficacy metrics that range between 0 and 1 (i.e. percent viability). Defaults to FALSE.
#' @param Average_Duplicate_Records A logic vector of length 1 indicating whether or not duplicated records (where a cell line has multiple records for being tested with a given drug at a given concentration) should be averaged. If TRUE, Efficacy values are averaged, and, if Calculate_Uncertainty is also TRUE, Efficacy_SE values are added in quadrature and divided by the number of duplicate records for that cell line/drug/concentration set.
#' @param Return_Bootstrap_Values A logic vector of length 1 indicating whether or not the function should return the Drug1 Efficacies and Drug2 Efficacies simulated in the semi-parametric bootstrap used to estimate the uncertainties of those values. If equal to TRUE, Calculate_Uncertainty must also equal TRUE.
#'
#'@details
#'Uncertainty estimates for values calculated by this function are generated using a semi-parametric bootstrap approach. This is performed in several steps.\enumerate{
#'\item Drug1 efficacies for each concentration are simulated by random sampling from normal distributions with means equal to the provided calculated efficacies and standard deviations equal to the provided efficacy standard errors.
#'\item Drug2 efficacies are simulated in the same fashion as Drug1 efficacies, except in cases when Drug1 equals Drug2. In such cases, it is assumed that the efficacy values for Drug1 and Drug2 are derived from the same dose-response curve, so each simulated efficacy for Drug2 is matched to the corresponding simulated efficacy from Drug1 using a standard normal deviate.
#'\item Efficacy predictions are made for the combination of Drug1 + Drug2 for each cell line and set of simulated efficacies using the assumptions of independent drug action.
#'\item Cell lines are randomly sampled with replacement for each simulation as many times as there are original cell lines. The simulated Drug1 monotherapy efficacies and Drug1+Drug2 combination efficacies are then sampled according to the sampled cell lines for each simulation.
#'\item Mean efficacies are calculated for the monotherapy and combination treatments for each simulation. If specified to do so, these values are then used to calculate simulated HRs and IDAcomboscores.
#'\item The simulated distributions of each efficacy metric are used to estimate uncertainties for those metrics.
#'}
#'
#'@return \itemize{
#'\item If Return_Bootstrap_Values = FALSE, this function returns a list with 4 elements: 1) Either a data frame with the calculated efficacy predictions, or, if an error occurred, a character vector of length one with the error message. 2) A character value with the name of Drug1 3) A character value with the name of Drug1 4) A character vector containing the names of the cell lines used to make the efficacy predictions.
#'\item If Return_Bootstrap_Values = TRUE & Calculate_Uncertainty = TRUE, this function returns a list with 6 elements: the first 4 elements are the same as when Return_Bootstrap_Values = FALSE and the fifth and sixth elements are numeric vectors of, respectively, the Drug1 and Drug2 viabilities simulated during the semi-parametric bootstrap used to estimate uncertainties.
#'}
#'
#' @examples
#' #Loading Package
#'   library(IDACombo)
#'
#' #Making fake monotherapy dataset
#'   CellLineNames <- rep(c("CL1", "CL2", "CL3", "CL4", "CL5", "CL6"), 4)
#'   DrugNames <- c(rep("D1", 12), rep("D2", 12))
#'   Concentrations <- c(rep(1, 6), rep(2, 6), rep("1.5", 6), rep("3", 6))
#'   Viability <- c(sample(seq(0.4,1,length.out = 10), 6, replace = TRUE),
#'                  sample(seq(0.2,0.8,length.out = 10), 6, replace = TRUE),
#'                  sample(seq(0.4,1,length.out = 10), 6, replace = TRUE),
#'                  sample(seq(0.2,0.6,length.out = 10), 6, replace = TRUE))
#'   Viability_SE <- Viability * sample(seq(0,0.1,length.out = 100), 24, replace = TRUE)
#'   Fake_Data <- data.frame(CellLineNames, DrugNames, Concentrations, Viability, Viability_SE)
#'
#' #Creating efficacy predictions for D1 + D2 without uncertainty calculations
#'   IDAPredict.2drug(Monotherapy_Data = Fake_Data,
#'                    Cell_Line_Name_Column = "CellLineNames",
#'                    Drug_Name_Column = "DrugNames",
#'                    Drug_Concentration_Column = "Concentrations",
#'                    Efficacy_Column = "Viability",
#'                    Drug1 = "D1",
#'                    Drug2 = "D2",
#'                    Calculate_Uncertainty = FALSE,
#'                    LowerEfficacyIsBetterDrugEffect = TRUE,
#'                    Efficacy_Metric_Name = "Viability",
#'                    Calculate_IDAcomboscore_And_Hazard_Ratio = TRUE,
#'                    Average_Duplicate_Records = FALSE)
#'
#' #Creating efficacy predictions for D1 + D2 with uncertainty calculations
#'   IDAPredict.2drug(Monotherapy_Data = Fake_Data,
#'                    Cell_Line_Name_Column = "CellLineNames",
#'                    Drug_Name_Column = "DrugNames",
#'                    Drug_Concentration_Column = "Concentrations",
#'                    Efficacy_Column = "Viability",
#'                    Drug1 = "D1",
#'                    Drug2 = "D2",
#'                    Calculate_Uncertainty = TRUE,
#'                    Efficacy_SE_Column = "Viability_SE",
#'                    LowerEfficacyIsBetterDrugEffect = TRUE,
#'                    Efficacy_Metric_Name = "Viability",
#'                    Calculate_IDAcomboscore_And_Hazard_Ratio = TRUE,
#'                    Average_Duplicate_Records = FALSE)
#'
#' #Creating efficacy predictions for D1 + D2 with uncertainty calculations
#' #and returning simulated values from semi-parametric bootstrap
#'   IDAPredict.2drug(Monotherapy_Data = Fake_Data,
#'                    Cell_Line_Name_Column = "CellLineNames",
#'                    Drug_Name_Column = "DrugNames",
#'                    Drug_Concentration_Column = "Concentrations",
#'                    Efficacy_Column = "Viability",
#'                    Drug1 = "D1",
#'                    Drug2 = "D2",
#'                    Calculate_Uncertainty = TRUE,
#'                    Efficacy_SE_Column = "Viability_SE",
#'                    LowerEfficacyIsBetterDrugEffect = TRUE,
#'                    Efficacy_Metric_Name = "Viability",
#'                    Calculate_IDAcomboscore_And_Hazard_Ratio = TRUE,
#'                    Average_Duplicate_Records = FALSE,
#'                    Return_Bootstrap_Values = TRUE)
#'
#' #Converting Viabilty to reduction in viability and redoing calculations
#' #Note the change in the LowerEfficacyIsBetterDrugEffect flag from TRUE to FALSE
#'   Reduction_in_Viability <- 1-Viability
#'   Reduction_in_Viability_SE <- Viability_SE
#'   Fake_Data <- data.frame(CellLineNames,
#'                           DrugNames,
#'                           Concentrations,
#'                           Reduction_in_Viability,
#'                           Reduction_in_Viability_SE)
#'   IDAPredict.2drug(Monotherapy_Data = Fake_Data,
#'                    Cell_Line_Name_Column = "CellLineNames",
#'                    Drug_Name_Column = "DrugNames",
#'                    Drug_Concentration_Column = "Concentrations",
#'                    Efficacy_Column = "Reduction_in_Viability",
#'                    Drug1 = "D1",
#'                    Drug2 = "D2",
#'                    Calculate_Uncertainty = TRUE,
#'                    LowerEfficacyIsBetterDrugEffect = FALSE,
#'                    Efficacy_SE_Column = "Reduction_in_Viability_SE",
#'                    Efficacy_Metric_Name = "Reduction_In_Viability",
#'                    Calculate_IDAcomboscore_And_Hazard_Ratio = TRUE,
#'                    Average_Duplicate_Records = FALSE)
#'
#' #Changing efficacy metric to percent growth (range -1 to 1)
#' #Note that calculating Hazard Ratios and IDAcomboscores is no longer valid, so
#' #Calculate_IDAcomboscore_And_Hazard_Ratio is set to FALSE.
#'   Percent_Growth <- c(sample(seq(0.4,1,length.out = 10), 6, replace = TRUE),
#'                       sample(seq(-0.4,0.2,length.out = 10), 6, replace = TRUE),
#'                       sample(seq(-0.2,0.3,length.out = 10), 6, replace = TRUE),
#'                       sample(seq(-1,0.2,length.out = 10), 6, replace = TRUE))
#'   Percent_Growth_SE <- abs(Percent_Growth * sample(seq(0,0.1,length.out = 100), 24, replace = TRUE))
#'   Fake_Data <- data.frame(CellLineNames,
#'                           DrugNames,
#'                           Concentrations,
#'                           Percent_Growth,
#'                           Percent_Growth_SE)
#'   IDAPredict.2drug(Monotherapy_Data = Fake_Data,
#'                    Cell_Line_Name_Column = "CellLineNames",
#'                    Drug_Name_Column = "DrugNames",
#'                    Drug_Concentration_Column = "Concentrations",
#'                    Efficacy_Column = "Percent_Growth",
#'                    Drug1 = "D1",
#'                    Drug2 = "D2",
#'                    Calculate_Uncertainty = TRUE,
#'                    LowerEfficacyIsBetterDrugEffect = TRUE,
#'                    Efficacy_SE_Column = "Percent_Growth_SE",
#'                    Calculate_IDAcomboscore_And_Hazard_Ratio = FALSE,
#'                    Efficacy_Metric_Name = "Percent_Growth",
#'                    Average_Duplicate_Records = FALSE)
#'
#' #Adding duplicate records for each cell line, and showing behavior with
#' #Average_Duplicate_Records = FALSE. Should produce warning messages that
#' #duplicates were found and removed.
#'   Percent_Growth <- c(sample(seq(0.4,1,length.out = 10), 6, replace = TRUE),
#'                       sample(seq(-0.4,0.2,length.out = 10), 6, replace = TRUE),
#'                       sample(seq(-0.2,0.3,length.out = 10), 6, replace = TRUE),
#'                       sample(seq(-1,0.2,length.out = 10), 6, replace = TRUE))
#'   Percent_Growth_SE <- abs(Percent_Growth * sample(seq(0,0.1,length.out = 100), 24, replace = TRUE))
#'   Fake_Data_to_add <- data.frame(CellLineNames,
#'                           DrugNames,
#'                           Concentrations,
#'                           Percent_Growth,
#'                           Percent_Growth_SE)
#'   Fake_Data <- rbind(Fake_Data, Fake_Data_to_add)
#'   IDAPredict.2drug(Monotherapy_Data = Fake_Data,
#'                    Cell_Line_Name_Column = "CellLineNames",
#'                    Drug_Name_Column = "DrugNames",
#'                    Drug_Concentration_Column = "Concentrations",
#'                    Efficacy_Column = "Percent_Growth",
#'                    Drug1 = "D1",
#'                    Drug2 = "D2",
#'                    Calculate_Uncertainty = TRUE,
#'                    LowerEfficacyIsBetterDrugEffect = TRUE,
#'                    Efficacy_SE_Column = "Percent_Growth_SE",
#'                    Calculate_IDAcomboscore_And_Hazard_Ratio = FALSE,
#'                    Efficacy_Metric_Name = "Percent_Growth",
#'                    Average_Duplicate_Records = FALSE)
#'
#' #Now setting to average duplicate values.
#'   Fake_Data <- rbind(Fake_Data, Fake_Data_to_add)
#'   IDAPredict.2drug(Monotherapy_Data = Fake_Data,
#'                    Cell_Line_Name_Column = "CellLineNames",
#'                    Drug_Name_Column = "DrugNames",
#'                    Drug_Concentration_Column = "Concentrations",
#'                    Efficacy_Column = "Percent_Growth",
#'                    Drug1 = "D1",
#'                    Drug2 = "D2",
#'                    Calculate_Uncertainty = TRUE,
#'                    LowerEfficacyIsBetterDrugEffect = TRUE,
#'                    Efficacy_SE_Column = "Percent_Growth_SE",
#'                    Calculate_IDAcomboscore_And_Hazard_Ratio = FALSE,
#'                    Efficacy_Metric_Name = "Percent_Growth",
#'                    Average_Duplicate_Records = TRUE)
#' @export
IDAPredict.2drug <- function(Monotherapy_Data, Cell_Line_Name_Column, Drug_Name_Column, Drug_Concentration_Column, Efficacy_Column, LowerEfficacyIsBetterDrugEffect, Efficacy_Metric_Name = "Efficacy", Drug1, Drug2, Calculate_Uncertainty = FALSE, Efficacy_SE_Column = NULL, n_Simulations = 1000, Calculate_IDAcomboscore_And_Hazard_Ratio = FALSE, Average_Duplicate_Records = FALSE, Return_Bootstrap_Values = FALSE){

  #Checking that all input variables are in correct format
    if(! is.data.frame(Monotherapy_Data)){
      stop("Monotherapy_Data is not a data frame.")
    }
    if(! is.vector(Cell_Line_Name_Column) | ! is.character(Cell_Line_Name_Column) | ! length(Cell_Line_Name_Column) == 1){
      stop("Cell_Line_Name_Column is not a character vector of length 1.")
    }
    if(! Cell_Line_Name_Column %in% colnames(Monotherapy_Data)){
      stop("Cell_Line_Name_Column does not match any column names in Monotherapy_Data.")
    }
    if(! is.vector(Drug_Name_Column) | ! is.character(Drug_Name_Column) | ! length(Drug_Name_Column) == 1){
      stop("Drug_Name_Column is not a character vector of length 1.")
    }
    if(! Drug_Name_Column %in% colnames(Monotherapy_Data)){
      stop("Drug_Name_Column does not match any column names in Monotherapy_Data.")
    }
    if(! is.vector(Drug_Concentration_Column) | ! is.character(Drug_Concentration_Column) | ! length(Drug_Concentration_Column) == 1){
      stop("Drug_Concentration_Column is not a character vector of length 1.")
    }
    if(! Drug_Concentration_Column %in% colnames(Monotherapy_Data)){
      stop("Drug_Concentration_Column does not match any column names in Monotherapy_Data.")
    }
    if(! is.vector(Efficacy_Column) | ! is.character(Efficacy_Column) | ! length(Efficacy_Column) == 1){
      stop("Efficacy_Column is not a character vector of length 1.")
    }
    if(! Efficacy_Column %in% colnames(Monotherapy_Data)){
      stop("Efficacy_Column does not match any column names in Monotherapy_Data.")
    }
    if(! is.vector(LowerEfficacyIsBetterDrugEffect) | ! is.logical(LowerEfficacyIsBetterDrugEffect) | ! length(LowerEfficacyIsBetterDrugEffect) == 1){
      stop("LowerEfficacyIsBetterDrugEffect is not a logical vector of length 1.")
    }
    if(! is.vector(Efficacy_Metric_Name) | ! is.character(Efficacy_Metric_Name) | ! length(Efficacy_Metric_Name) == 1){
      stop("Efficacy_Metric_Name is not a character vector of length 1.")
    }
    if(! is.vector(Drug1) | ! is.character(Drug1) | ! length(Drug1) == 1){
      stop("Drug1 is not a character vector with length 1.")
    }
    if(! Drug1 %in% Monotherapy_Data[,Drug_Name_Column]){
      stop(paste0("No ", Drug1, " data found in Monotherapy_Data."))
    }
    if(! is.vector(Drug2) | ! is.character(Drug2) | ! length(Drug2) == 1){
      stop("Drug2 is not a character vector with length 1.")
    }
    if(! Drug2 %in% Monotherapy_Data[,Drug_Name_Column]){
      stop(paste0("No ", Drug2, " data found in Monotherapy_Data."))
    }
    if(! is.vector(Calculate_Uncertainty) | ! is.logical(Calculate_Uncertainty) | ! length(Calculate_Uncertainty) == 1){
      stop("Calculate_Uncertainty is not a logical vector of length 1.")
    }
    if(Calculate_Uncertainty == TRUE){
      if(is.null(Efficacy_SE_Column)){
        stop("Calculate_Uncertainty is TRUE but Efficacy_SE_Column is not specified. Please either set Calculate_Uncertainty to FALSE or specify Efficacy_SE_Column.")
      }
      if(! is.vector(Efficacy_SE_Column) | ! is.character(Efficacy_SE_Column) | ! length(Efficacy_SE_Column) == 1){
        stop("Calculate_Uncertainty is TRUE, but Efficacy_SE_Column is not a character vector of length 1.")
      }
      if(! Efficacy_SE_Column %in% colnames(Monotherapy_Data)){
        stop("Calculate_Uncertainty is TRUE, but Efficacy_SE_Column does not match any column names in Monotherapy_Data.")
      }
      if(! is.vector(n_Simulations) | ! is.numeric(n_Simulations) | ! length(n_Simulations) == 1){
        stop("Calculate_Uncertainty is TRUE, but n_Simulations is not a numeric vector of length 1.")
      }
      if(! n_Simulations%%1==0 | ! n_Simulations >= 40){
        stop("Calculate_Uncertainty is TRUE, but n_Simulations is not a positive integer >= 40.")
      }
    }
    if(! is.vector(Return_Bootstrap_Values) | ! is.logical(Return_Bootstrap_Values) | ! length(Return_Bootstrap_Values) == 1){
      stop("Return_Bootstrap_Values is not a logical vector of length 1.")
    }
    if(Return_Bootstrap_Values == TRUE & ! Calculate_Uncertainty == TRUE){
      stop("Return_Bootstrap_Values is TRUE, but Calculate_Uncertainty is not TRUE.")
    }
    if(! is.vector(Calculate_IDAcomboscore_And_Hazard_Ratio) | ! is.logical(Calculate_IDAcomboscore_And_Hazard_Ratio) | ! length(Calculate_IDAcomboscore_And_Hazard_Ratio) == 1){
      stop("Calculate_IDAcomboscore_And_Hazard_Ratio is not a logical vector of length 1.")
    }
    if(! is.vector(Average_Duplicate_Records) | ! is.logical(Average_Duplicate_Records) | ! length(Average_Duplicate_Records) == 1){
      stop("Average_Duplicate_Records is not a logical vector of length 1.")
    }


  #Organizing data into standard format based on column names provided for each desired set of information
  #Also subsetting to only include data pertaining to Drug1 and Drug2
    if(Calculate_Uncertainty == TRUE){
      Data <- Monotherapy_Data[Monotherapy_Data[,Drug_Name_Column] %in% c(Drug1, Drug2),c(Cell_Line_Name_Column, Drug_Name_Column, Drug_Concentration_Column, Efficacy_Column, Efficacy_SE_Column)]
      colnames(Data) <- c("CellLine", "Drug", "Conc", "Efficacy", "Efficacy_SE")
    } else {
      Data <- Monotherapy_Data[Monotherapy_Data[,Drug_Name_Column] %in% c(Drug1, Drug2),c(Cell_Line_Name_Column, Drug_Name_Column, Drug_Concentration_Column, Efficacy_Column)]
      colnames(Data) <- c("CellLine", "Drug", "Conc", "Efficacy")
    }
    rm(Monotherapy_Data)

  #Making sure all columns are in correct formats
    Data$CellLine <- as.character(Data$CellLine)
    Data$Drug <- as.character(Data$Drug)
    Data$Conc <- as.character(Data$Conc)
    Data$Efficacy <- as.numeric(as.character(Data$Efficacy))
    if(Calculate_Uncertainty == TRUE){
      Data$Efficacy_SE <- as.numeric(as.character(Data$Efficacy_SE))
    }

  #Removing rows that are missing information
  #Note: missing SE information is ignored if Calculate_Uncertainty == FALSE
    Data <- Data[complete.cases(Data),]


  #Subsetting into Drug1 and Drug2 data
    Drug1Data <- Data[Data$Drug %in% Drug1,]
    Drug2Data <- Data[Data$Drug %in% Drug2,]
    rm(Data)

  #Finding cell line overlap between all drugs
    Usable_CellLines <- sort(unique(Drug1Data$CellLine[Drug1Data$CellLine %in% Drug2Data$CellLine]))

  #Checking if at least 2 cell lines overlap for all drugs. If not, exiting with no
  #predictions and a warning.
    if(! length(Usable_CellLines) >= 2){
      #Returning NA predictions with warning due to too few cell lines.
      warning(paste0("<2 overlapping cell lines available for combination of ", Drug1, " + ", Drug2, "."))
      if(Return_Bootstrap_Values == FALSE){
        Return_Object <- list("Less than 2 overlapping cell lines available.", Drug1, Drug2, Usable_CellLines)
        names(Return_Object) <- c("Efficacy_Predictions", "Drug1", "Drug2", "Cell_Lines_Used")
        return(Return_Object)
      } else if(Return_Bootstrap_Values == TRUE){
        Return_Object <- list("Less than 2 overlapping cell lines available.", Drug1, Drug2, Usable_CellLines, NULL, NULL)
        names(Return_Object) <- c("Efficacy_Predictions", "Drug1", "Drug2", "Cell_Lines_Used", "Bootstrap_Drug1_Efficacies", "Bootstrap_Drug2_Efficacies")
        return(Return_Object)
      }
    }

  #Subsetting drug data to only include overlapping cell lines and ordering same
  #for each drug.
    Drug1Data <- Drug1Data[Drug1Data$CellLine %in% Usable_CellLines,]
    Drug1Data <- Drug1Data[order(Drug1Data$Conc, Drug1Data$CellLine),]
    Drug2Data <- Drug2Data[Drug2Data$CellLine %in% Usable_CellLines,]
    Drug2Data <- Drug2Data[order(Drug2Data$Conc, Drug2Data$CellLine),]

  #Checking that cell lines aren't duplicated in each drug dataset for a given drug dose.
  #If Average_Duplicate_Records == FALSE, removing cell line duplicates with warning if duplicates are found.
  #If Average_Duplicate_Records == TRUE, averaging duplicate records without warning.
    D1_CL_Conc <- paste(Drug1Data$CellLine, Drug1Data$Conc, sep = "_")
    D1Dups <- D1_CL_Conc[duplicated(D1_CL_Conc)]
    if(length(D1Dups) > 0 & Average_Duplicate_Records == FALSE){
      warning(paste0("Duplicated information found for the following cell lines and ", Drug1, " concentrations. Average_Duplicate_Records = FALSE so duplicates removed: ", paste(D1Dups, collapse = ", ")))
      Drug1Data <- Drug1Data[! duplicated(D1_CL_Conc),]
    } else if(length(D1Dups) > 0 & Average_Duplicate_Records == TRUE){
      if(Calculate_Uncertainty == FALSE){
        #Simply averaging efficacy values
          d1.colnames <- colnames(Drug1Data)
          Drug1Data <- aggregate(Drug1Data$Efficacy, by = list(Drug1Data$CellLine, Drug1Data$Drug, Drug1Data$Conc), FUN = mean)
          colnames(Drug1Data) <- d1.colnames
      } else if(Calculate_Uncertainty == TRUE){
        #Averaging efficacy
          d1.Efficacy_Average <- aggregate(Drug1Data$Efficacy, by = list(Drug1Data$CellLine, Drug1Data$Drug, Drug1Data$Conc), FUN = mean)
          colnames(d1.Efficacy_Average) <- c("CellLine", "Drug", "Conc", "Efficacy")
        #Calculating uncertainty in averaged efficacy by adding efficacy uncertainties in quadrature and dividing by number of values used in average
          d1.Efficacy_Average_Uncertainties <- aggregate(Drug1Data$Efficacy_SE, by = list(Drug1Data$CellLine, Drug1Data$Drug, Drug1Data$Conc), FUN = function(x){sqrt(sum(x^2))/length(x)})
          colnames(d1.Efficacy_Average_Uncertainties) <- c("CellLine", "Drug", "Conc", "Efficacy_SE")
        #Combining results
          Drug1Data <- merge(d1.Efficacy_Average, d1.Efficacy_Average_Uncertainties)
      }
    }
    D2_CL_Conc <- paste(Drug2Data$CellLine, Drug2Data$Conc, sep = "_")
    D2Dups <- D2_CL_Conc[duplicated(D2_CL_Conc)]
    if(length(D2Dups) > 0 & Average_Duplicate_Records == FALSE){
      warning(paste0("Duplicated information found for the following cell lines and ", Drug2, " concentrations. Average_Duplicate_Records = FALSE so duplicates removed: ", paste(D2Dups, collapse = ", ")))
      Drug2Data <- Drug2Data[! duplicated(D2_CL_Conc),]
    } else if(length(D2Dups) > 0 & Average_Duplicate_Records == TRUE){
      if(Calculate_Uncertainty == FALSE){
        #Simply averaging efficacy values
          d2.colnames <- colnames(Drug2Data)
          Drug2Data <- aggregate(Drug2Data$Efficacy, by = list(Drug2Data$CellLine, Drug2Data$Drug, Drug2Data$Conc), FUN = mean)
          colnames(Drug2Data) <- d2.colnames
      } else if(Calculate_Uncertainty == TRUE){
        #Averaging efficacy
          d2.Efficacy_Average <- aggregate(Drug2Data$Efficacy, by = list(Drug2Data$CellLine, Drug2Data$Drug, Drug2Data$Conc), FUN = mean)
          colnames(d2.Efficacy_Average) <- c("CellLine", "Drug", "Conc", "Efficacy")
        #Calculating uncertainty in averaged efficacy by adding efficacy uncertainties in quadrature and dividing by number of values used in average
          d2.Efficacy_Average_Uncertainties <- aggregate(Drug2Data$Efficacy_SE, by = list(Drug2Data$CellLine, Drug2Data$Drug, Drug2Data$Conc), FUN = function(x){sqrt(sum(x^2))/length(x)})
          colnames(d2.Efficacy_Average_Uncertainties) <- c("CellLine", "Drug", "Conc", "Efficacy_SE")
        #Combining results
          Drug2Data <- merge(d2.Efficacy_Average, d2.Efficacy_Average_Uncertainties)
      }
    }

  #Checking that all drug concentrations for each drug are available for each cell line
  #Omitting cell lines that are missing concentrations for 1 or more drugs
    AllData <- list(Drug1Data, Drug2Data)
    CLs_per_dose <- lapply(AllData, function(x){as.data.frame.table(table(x$Conc), stringsAsFactors = FALSE)[,2]})
    if(! all(CLs_per_dose[[1]] == length(Usable_CellLines)) | ! all(CLs_per_dose[[2]] == length(Usable_CellLines))){
      #Not all cell lines have all doses for all drugs. Printing warning and removing cell lines with incomplete information.
        for(i in 1:length(AllData)){
          n_doses <- length(unique(AllData[[i]]$Conc))
          CL_dose_count <- as.data.frame.table(table(AllData[[i]]$CellLine), stringsAsFactors = FALSE)
          if(i == 1){
            CLs_missing_doses <- CL_dose_count$Var1[CL_dose_count$Freq < n_doses]
          } else {
            CLs_missing_doses <- unique(c(CLs_missing_doses, CL_dose_count$Var1[CL_dose_count$Freq < n_doses]))
          }
        }
        warning(paste0("The following cell lines are missing efficacy information for one or more drug concentrations in the combination of ", Drug1, " + ", Drug2, ". These cell lines have been omitted from the analysis: ", paste(CLs_missing_doses, collapse = ", ")))
      #Re-subsetting drug data to only include overlapping cell lines and ordering same
      #for each drug with cell lines that had missing information removed.
        Usable_CellLines <- Usable_CellLines[! Usable_CellLines %in% CLs_missing_doses]
        Drug1Data <- Drug1Data[Drug1Data$CellLine %in% Usable_CellLines,]
        Drug1Data <- Drug1Data[order(Drug1Data$Conc, Drug1Data$CellLine),]
        Drug2Data <- Drug2Data[Drug2Data$CellLine %in% Usable_CellLines,]
        Drug2Data <- Drug2Data[order(Drug2Data$Conc, Drug2Data$CellLine),]
        rm(AllData, CLs_per_dose, n_doses, CLs_missing_doses, CL_dose_count)
    } else {
      rm(AllData, CLs_per_dose)
    }

  #Checking if at least 2 cell lines overlap for all drugs. If not, exiting with no
  #predictions and a warning.
    if(! length(Usable_CellLines) >= 2){
      #Returning NA predictions with warning due to too few cell lines.
      warning(paste0("<2 overlapping cell lines available for combination of ", Drug1, " + ", Drug2, "."))
      if(Return_Bootstrap_Values == FALSE){
        Return_Object <- list("Less than 2 overlapping cell lines available.", Drug1, Drug2, Usable_CellLines)
        names(Return_Object) <- c("Efficacy_Predictions", "Drug1", "Drug2", "Cell_Lines_Used")
        return(Return_Object)
      } else if(Return_Bootstrap_Values == TRUE){
        Return_Object <- list("Less than 2 overlapping cell lines available.", Drug1, Drug2, Usable_CellLines, NULL, NULL)
        names(Return_Object) <- c("Efficacy_Predictions", "Drug1", "Drug2", "Cell_Lines_Used", "Bootstrap_Drug1_Efficacies", "Bootstrap_Drug2_Efficacies")
        return(Return_Object)
      }
    }

  #Dividing data by drug concentration
    Drug1Data <- split(Drug1Data, as.factor(Drug1Data$Conc))
    Drug2Data <- split(Drug2Data, as.factor(Drug2Data$Conc))

  #Identifying all concentration comparisons that need to be made.
    D1_concentrations <- names(Drug1Data)
    D2_concentrations <- names(Drug2Data)
    Dose_Comparisons <- expand.grid(D1_concentrations, D2_concentrations, stringsAsFactors = FALSE)
    colnames(Dose_Comparisons) <- c("Drug1Dose", "Drug2Dose")

  #Adding extra columns to Dose_Comparisons to store predicted efficacy results
    Dose_Comparisons$Mean_Drug1_Efficacy <- NA
    Dose_Comparisons$Mean_Drug1_Efficacy_SE <- NA
    Dose_Comparisons$`Mean_Drug1_Efficacy_95%_Confidence_Interval` <- NA
    Dose_Comparisons$Mean_Drug2_Efficacy <- NA
    Dose_Comparisons$Mean_Drug2_Efficacy_SE <- NA
    Dose_Comparisons$`Mean_Drug2_Efficacy_95%_Confidence_Interval` <- NA
    Dose_Comparisons$Mean_Combo_Efficacy <- NA
    Dose_Comparisons$Mean_Combo_Efficacy_SE <- NA
    Dose_Comparisons$`Mean_Combo_Efficacy_95%_Confidence_Interval` <- NA
    Dose_Comparisons$HR_vs_Drug1 <- NA
    Dose_Comparisons$HR_vs_Drug1_SE <- NA
    Dose_Comparisons$`HR_vs_Drug1_95%_Confidence_Interval` <- NA
    Dose_Comparisons$`p_HR_vs_Drug1>=1` <- NA
    # Dose_Comparisons$`p_HR_vs_Drug1=1` <- NA
    Dose_Comparisons$HR_vs_Drug2 <- NA
    Dose_Comparisons$HR_vs_Drug2_SE <- NA
    Dose_Comparisons$`HR_vs_Drug2_95%_Confidence_Interval` <- NA
    Dose_Comparisons$`p_HR_vs_Drug2>=1` <- NA
    # Dose_Comparisons$`p_HR_vs_Drug2=1` <- NA
    Dose_Comparisons$IDA_Comboscore <- NA
    Dose_Comparisons$IDA_Comboscore_SE <- NA
    Dose_Comparisons$`IDA_Comboscore_95%_Confidence_Interval` <- NA
    Dose_Comparisons$`p_IDA_Comboscore<=0` <- NA

  #Performing combination efficacy predictions for cases where lower efficacy values
  #indicate a more effective drug effect (i.e. efficacy = percent viability, percent growth, etc.)
    if(LowerEfficacyIsBetterDrugEffect == TRUE){
      #Looping through all dose comparisons and predicting combination efficacies
        for(i in 1:nrow(Dose_Comparisons)){
          #Identifying correct dose data for each drug for this comparison.
            D1_Data <- Drug1Data[[which(names(Drug1Data) == Dose_Comparisons$Drug1Dose[i])]]
            D2_Data <- Drug2Data[[which(names(Drug2Data) == Dose_Comparisons$Drug2Dose[i])]]
          #Calculating expected combination efficacy for each cell line
          #using Independent Drug Action
            Combo_Efficacy <- pmin(D1_Data$Efficacy, D2_Data$Efficacy)
          #Calculating average efficacy across all cell lines
            Dose_Comparisons$Mean_Drug1_Efficacy[i] <- mean(D1_Data$Efficacy)
            Dose_Comparisons$Mean_Drug2_Efficacy[i] <- mean(D2_Data$Efficacy)
            Dose_Comparisons$Mean_Combo_Efficacy[i] <- mean(Combo_Efficacy)
            if(Calculate_IDAcomboscore_And_Hazard_Ratio == TRUE){
              Dose_Comparisons$HR_vs_Drug1[i] <- Dose_Comparisons$Mean_Combo_Efficacy[i] / Dose_Comparisons$Mean_Drug1_Efficacy[i]
              Dose_Comparisons$HR_vs_Drug2[i] <- Dose_Comparisons$Mean_Combo_Efficacy[i] / Dose_Comparisons$Mean_Drug2_Efficacy[i]
              delta_viability <- min(Dose_Comparisons$Mean_Drug1_Efficacy[i], Dose_Comparisons$Mean_Drug2_Efficacy[i]) - Dose_Comparisons$Mean_Combo_Efficacy[i]
              HR_C_over_Mbest <- max(Dose_Comparisons$HR_vs_Drug1[i], Dose_Comparisons$HR_vs_Drug2[i])
              Dose_Comparisons$IDA_Comboscore[i] <- delta_viability - delta_viability * HR_C_over_Mbest
            }
          }
        rm(D1_Data, D2_Data, Combo_Efficacy)
        if(Calculate_IDAcomboscore_And_Hazard_Ratio == TRUE){
          rm(delta_viability, HR_C_over_Mbest)
        }

      #If Calculate_Uncertainty == TRUE, doing Monte Carlo simulation to estimate
      #uncertainties in output parameters based on uncertainties in monotherapy efficacies.
        if(Calculate_Uncertainty == TRUE){
          #Looping through each drug for each treatment and simulating efficacies based on
          #measured efficacies and SE's
            D1_MC_Efficacies <- as.list(NULL)
            for(i in 1:length(Drug1Data)){
              if(i == 1){
                D1_MC_Efficacies[[i]] <- apply(Drug1Data[[i]][,c("Efficacy", "Efficacy_SE")], 1, function(x){rnorm(n = n_Simulations, mean = x[1], sd = x[2])})
                colnames(D1_MC_Efficacies[[i]]) <- Drug1Data[[i]]$CellLine
                Measured_D1_Data <- Drug1Data[[names(Drug1Data)[i]]][,c("CellLine", "Efficacy", "Efficacy_SE")]
                Measured_D1_Data <- Measured_D1_Data[match(colnames(D1_MC_Efficacies[[1]]), Measured_D1_Data$CellLine),]
                Measured_D1_Efficacies <- matrix(Measured_D1_Data$Efficacy, ncol = length(Measured_D1_Data$Efficacy), nrow = n_Simulations, byrow = TRUE)
                Measured_D1_Efficacy_SEs <- matrix(Measured_D1_Data$Efficacy_SE, ncol = length(Measured_D1_Data$Efficacy_SE), nrow = n_Simulations, byrow = TRUE)
                SEs_deviated <- (D1_MC_Efficacies[[1]] - Measured_D1_Efficacies) / Measured_D1_Efficacy_SEs
              } else {
                Measured_D1_Data <- Drug1Data[[names(Drug1Data)[i]]][,c("CellLine", "Efficacy", "Efficacy_SE")]
                Measured_D1_Data <- Measured_D1_Data[match(colnames(SEs_deviated), Measured_D1_Data$CellLine),]
                Measured_D1_Efficacies <- matrix(Measured_D1_Data$Efficacy, ncol = length(Measured_D1_Data$Efficacy), nrow = n_Simulations, byrow = TRUE, dimnames = list(NULL, Measured_D1_Data$CellLine))
                Measured_D1_Efficacy_SEs <- matrix(Measured_D1_Data$Efficacy_SE, ncol = length(Measured_D1_Data$Efficacy_SE), nrow = n_Simulations, byrow = TRUE, dimnames = list(NULL, Measured_D1_Data$CellLine))
                D1_MC_Efficacies[[i]] <- Measured_D1_Efficacies + (SEs_deviated * Measured_D1_Efficacy_SEs)
              }
            }
            names(D1_MC_Efficacies) <- names(Drug1Data)
            D2_MC_Efficacies <- as.list(NULL)
            for(i in 1:length(Drug2Data)){
              #Checking if drug 2 is the same as drug 1.
              #If so, not randomly sampling for D2 viabilities. Instead, calculating how many SE's from
              #the measured value each simulated D1 viability fell, and matching that distance in
              #the simulated D2 viabilities. This is done because D1 and D2 viabilities are not independent in
              #this case--they will have come from the same dose-response curve.
                if(Drug1 == Drug2){
                  #Calculating simulated D2 viabilities based on the SE deviations from the simulated D1 viabilities
                    Measured_D2_Data <- Drug2Data[[names(Drug2Data)[i]]][,c("CellLine", "Efficacy", "Efficacy_SE")]
                    Measured_D2_Data <- Measured_D2_Data[match(colnames(SEs_deviated), Measured_D2_Data$CellLine),]
                    Measured_D2_Efficacies <- matrix(Measured_D2_Data$Efficacy, ncol = length(Measured_D2_Data$Efficacy), nrow = n_Simulations, byrow = TRUE, dimnames = list(NULL, Measured_D2_Data$CellLine))
                    Measured_D2_Efficacy_SEs <- matrix(Measured_D2_Data$Efficacy_SE, ncol = length(Measured_D2_Data$Efficacy_SE), nrow = n_Simulations, byrow = TRUE, dimnames = list(NULL, Measured_D2_Data$CellLine))
                    D2_MC_Efficacies[[i]] <- Measured_D2_Efficacies + (SEs_deviated * Measured_D2_Efficacy_SEs)
                } else {
                  #If Drug1 does not equal Drug2, randomly sampling
                    D2_MC_Efficacies[[i]] <- apply(Drug2Data[[i]][,c("Efficacy", "Efficacy_SE")], 1, function(x){rnorm(n = n_Simulations, mean = x[1], sd = x[2])})
                }
            }
            names(D2_MC_Efficacies) <- names(Drug2Data)
            #Cleaning up
              rm(Measured_D1_Data, Measured_D1_Efficacies, Measured_D1_Efficacy_SEs, SEs_deviated)
              if(exists("Measured_D2_Data")){
                rm(Measured_D2_Data, Measured_D2_Efficacies, Measured_D2_Efficacy_SEs)
              }
          #Resampling cell lines for each simulation with replacement to account for the effect of random variation in cell line selection
            Selected_CCLs_Per_Sim <- t(apply(D1_MC_Efficacies[[1]], 1, function(x){return(sample(1:length(x), replace = TRUE))}))
            for(i in 1:length(D1_MC_Efficacies)){
              D1_MC_Efficacies[[i]] <- t(mapply(function(x,y){return(x[y])}, x = split(D1_MC_Efficacies[[i]], row(D1_MC_Efficacies[[i]])), y = split(Selected_CCLs_Per_Sim, row(Selected_CCLs_Per_Sim))))
            }
            for(i in 1:length(D2_MC_Efficacies)){
              D2_MC_Efficacies[[i]] <- t(mapply(function(x,y){return(x[y])}, x = split(D2_MC_Efficacies[[i]], row(D2_MC_Efficacies[[i]])), y = split(Selected_CCLs_Per_Sim, row(Selected_CCLs_Per_Sim))))
            }
          #Looping through all dose comparisons and calculating uncertainties in output values
            for(i in 1:nrow(Dose_Comparisons)){
              #Identifying correct dose data for each drug for this comparison.
                D1_Data <- D1_MC_Efficacies[[which(names(D1_MC_Efficacies) == Dose_Comparisons$Drug1Dose[i])]]
                D2_Data <- D2_MC_Efficacies[[which(names(D2_MC_Efficacies) == Dose_Comparisons$Drug2Dose[i])]]
              #Calculating expected combination efficacy for each cell line
              #using Independent Drug Action
                Combo_Efficacy <- pmin(D1_Data, D2_Data)
              #Calculating average efficacy across sampled cell lines
                D1_efficacies <- rowMeans(D1_Data)
                D2_efficacies <- rowMeans(D2_Data)
                Combo_efficacies <- rowMeans(Combo_Efficacy)
              #Calculating standard errors
                Dose_Comparisons$Mean_Drug1_Efficacy_SE[i] <- sd(D1_efficacies)
                Dose_Comparisons$Mean_Drug2_Efficacy_SE[i] <- sd(D2_efficacies)
                Dose_Comparisons$Mean_Combo_Efficacy_SE[i] <- sd(Combo_efficacies)
              #Calculating 95% confidence intervals for mean treatment efficacies
                index_2.5 <- floor(0.025*n_Simulations)
                index_97.5 <- ceiling(0.975*n_Simulations)
                Dose_Comparisons$`Mean_Drug1_Efficacy_95%_Confidence_Interval`[i] <- paste(sort(D1_efficacies, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                Dose_Comparisons$`Mean_Drug2_Efficacy_95%_Confidence_Interval`[i] <- paste(sort(D2_efficacies, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                Dose_Comparisons$`Mean_Combo_Efficacy_95%_Confidence_Interval`[i] <- paste(sort(Combo_efficacies, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
              #Calculating simulated IDAcomboscores and HRs
                if(Calculate_IDAcomboscore_And_Hazard_Ratio == TRUE){
                  #HRs
                    HRs_vs_D1 <- Combo_efficacies / D1_efficacies
                    HRs_vs_D2 <- Combo_efficacies / D2_efficacies
                  #HR standard errors
                    Dose_Comparisons$HR_vs_Drug1_SE[i] <- sd(HRs_vs_D1)
                    Dose_Comparisons$HR_vs_Drug2_SE[i] <- sd(HRs_vs_D2)
                  #HR 95% confidence intervals
                    Dose_Comparisons$`HR_vs_Drug1_95%_Confidence_Interval`[i] <- paste(sort(HRs_vs_D1, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                    Dose_Comparisons$`HR_vs_Drug2_95%_Confidence_Interval`[i] <- paste(sort(HRs_vs_D2, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                  #Calculating HR p values
                    #Drug 1
                      # n_equal_to_1 <- sum(HRs_vs_D1 == 1)
                      # n_less_than_1 <- sum(HRs_vs_D1 < 1)
                      # n_greater_than_1 <- sum(HRs_vs_D1 > 1)
                      # lower_tail <- (n_less_than_1 + 0.5*n_equal_to_1)/length(HRs_vs_D1)
                      # upper_tail <- (n_greater_than_1 + 0.5*n_equal_to_1)/length(HRs_vs_D1)
                      #Calculating one-sided p-value with null hypothesis that HR >= 1
                        Dose_Comparisons$`p_HR_vs_Drug1>=1`[i] <- sum(HRs_vs_D1 >= 1) / length(HRs_vs_D1)
                        #If p-value is 0, setting as p < minimum p value that can be estimated using this many simulations
                          if(Dose_Comparisons$`p_HR_vs_Drug1>=1`[i] == 0){Dose_Comparisons$`p_HR_vs_Drug1>=1`[i] <- paste0("<", 1/n_Simulations)}
                      # #Calculating two-sided p-value with null hypothesis that HR = 1
                      #   Dose_Comparisons$`p_HR_vs_Drug1=1`[i] <- min(lower_tail, upper_tail)*2
                      #   #If p-value is 0, setting as p < minimum p value that can be estimated using this many simulations
                      #     if(Dose_Comparisons$`p_HR_vs_Drug1=1`[i] == 0){Dose_Comparisons$`p_HR_vs_Drug1=1`[i] <- paste0("<", 1/n_Simulations)}
                    #Drug 2
                      # n_equal_to_1 <- sum(HRs_vs_D2 == 1)
                      # n_less_than_1 <- sum(HRs_vs_D2 < 1)
                      # n_greater_than_1 <- sum(HRs_vs_D2 > 1)
                      # lower_tail <- (n_less_than_1 + 0.5*n_equal_to_1)/length(HRs_vs_D2)
                      # upper_tail <- (n_greater_than_1 + 0.5*n_equal_to_1)/length(HRs_vs_D2)
                      #Calculating one-sided p-value with null hypothesis that HR >= 1
                        Dose_Comparisons$`p_HR_vs_Drug2>=1`[i] <- sum(HRs_vs_D2 >= 1) / length(HRs_vs_D2)
                        #If p-value is 0, setting as p < minimum p value that can be estimated using this many simulations
                          if(Dose_Comparisons$`p_HR_vs_Drug2>=1`[i] == 0){Dose_Comparisons$`p_HR_vs_Drug2>=1`[i] <- paste0("<", 1/n_Simulations)}
                      # #Calculating two-sided p-value with null hypothesis that HR = 1
                      #   Dose_Comparisons$`p_HR_vs_Drug2=1`[i] <- min(lower_tail, upper_tail)*2
                      #   #If p-value is 0, setting as p < minimum p value that can be estimated using this many simulations
                      #     if(Dose_Comparisons$`p_HR_vs_Drug2=1`[i] == 0){Dose_Comparisons$`p_HR_vs_Drug2=1`[i] <- paste0("<", 1/n_Simulations)}
                  #IDAcomboscores
                    delta_viabilities <- pmin(D1_efficacies, D2_efficacies) - Combo_efficacies
                    HR_C_over_Mbests <- pmax(HRs_vs_D1, HRs_vs_D2)
                    MC_IDAcomboscores <- delta_viabilities - delta_viabilities * HR_C_over_Mbests
                  #IDAcomboscore standard error
                    Dose_Comparisons$IDA_Comboscore_SE[i] <- sd(MC_IDAcomboscores)
                  #IDAcomboscore 95% confidence interval
                    Dose_Comparisons$`IDA_Comboscore_95%_Confidence_Interval`[i] <- paste(sort(MC_IDAcomboscores, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                  #Calculating IDAcomboscore p-value
                    Dose_Comparisons$`p_IDA_Comboscore<=0`[i] <- sum(MC_IDAcomboscores <= 0) / length(MC_IDAcomboscores)
                    #If p-value is 0, setting as p < minimum p value that can be estimated using this many simulations
                      if(Dose_Comparisons$`p_IDA_Comboscore<=0`[i] == 0){Dose_Comparisons$`p_IDA_Comboscore<=0`[i] <- paste0("<", 1/n_Simulations)}
                }
            }
          #Cleaning up
            rm(D1_Data, D2_Data, Combo_Efficacy, index_2.5, index_97.5, D1_efficacies, D2_efficacies, Combo_efficacies)
            if(Calculate_IDAcomboscore_And_Hazard_Ratio == TRUE){
              # rm(delta_viabilities, HR_C_over_Mbests, n_equal_to_1, n_less_than_1, n_greater_than_1, lower_tail, upper_tail)
              rm(delta_viabilities, HR_C_over_Mbests, MC_IDAcomboscores, HRs_vs_D1, HRs_vs_D2)
            }
        }
    }

  #Performing combination efficacy predictions for cases where lower efficacy values
  #indicate a less effective drug effect (i.e. efficacy = percent cell death, etc.)
    if(LowerEfficacyIsBetterDrugEffect == FALSE){
      #Looping through all dose comparisons and predicting combination efficacies
        for(i in 1:nrow(Dose_Comparisons)){
          #Identifying correct dose data for each drug for this comparison.
            D1_Data <- Drug1Data[[which(names(Drug1Data) == Dose_Comparisons$Drug1Dose[i])]]
            D2_Data <- Drug2Data[[which(names(Drug2Data) == Dose_Comparisons$Drug2Dose[i])]]
          #Calculating expected combination efficacy for each cell line
          #using Independent Drug Action
            Combo_Efficacy <- pmax(D1_Data$Efficacy, D2_Data$Efficacy)
          #Calculating average efficacy across all cell lines
            Dose_Comparisons$Mean_Drug1_Efficacy[i] <- mean(D1_Data$Efficacy)
            Dose_Comparisons$Mean_Drug2_Efficacy[i] <- mean(D2_Data$Efficacy)
            Dose_Comparisons$Mean_Combo_Efficacy[i] <- mean(Combo_Efficacy)
            if(Calculate_IDAcomboscore_And_Hazard_Ratio == TRUE){
              Dose_Comparisons$HR_vs_Drug1[i] <- (1 - Dose_Comparisons$Mean_Combo_Efficacy[i]) / (1 - Dose_Comparisons$Mean_Drug1_Efficacy[i])
              Dose_Comparisons$HR_vs_Drug2[i] <- (1 - Dose_Comparisons$Mean_Combo_Efficacy[i]) / (1 - Dose_Comparisons$Mean_Drug2_Efficacy[i])
              delta_hazard <- min((1 - Dose_Comparisons$Mean_Drug1_Efficacy[i]), (1 - Dose_Comparisons$Mean_Drug2_Efficacy[i])) - (1 - Dose_Comparisons$Mean_Combo_Efficacy[i])
              HR_C_over_Mbest <- max(Dose_Comparisons$HR_vs_Drug1[i], Dose_Comparisons$HR_vs_Drug2[i])
              Dose_Comparisons$IDA_Comboscore[i] <- delta_hazard - delta_hazard * HR_C_over_Mbest
            }
          }
        rm(D1_Data, D2_Data, Combo_Efficacy)
        if(Calculate_IDAcomboscore_And_Hazard_Ratio == TRUE){
          rm(delta_hazard, HR_C_over_Mbest)
        }

      #If Calculate_Uncertainty == TRUE, doing Monte Carlo simulation to estimate
      #uncertainties in output parameters based on uncertainties in monotherapy efficacies.
        if(Calculate_Uncertainty == TRUE){
          #Looping through each drug for each treatment and simulating efficacies based on
          #measured efficacies and SE's
            D1_MC_Efficacies <- as.list(NULL)
            for(i in 1:length(Drug1Data)){
              if(i == 1){
                D1_MC_Efficacies[[i]] <- apply(Drug1Data[[i]][,c("Efficacy", "Efficacy_SE")], 1, function(x){rnorm(n = n_Simulations, mean = x[1], sd = x[2])})
                colnames(D1_MC_Efficacies[[i]]) <- Drug1Data[[i]]$CellLine
                Measured_D1_Data <- Drug1Data[[names(Drug1Data)[i]]][,c("CellLine", "Efficacy", "Efficacy_SE")]
                Measured_D1_Data <- Measured_D1_Data[match(colnames(D1_MC_Efficacies[[1]]), Measured_D1_Data$CellLine),]
                Measured_D1_Efficacies <- matrix(Measured_D1_Data$Efficacy, ncol = length(Measured_D1_Data$Efficacy), nrow = n_Simulations, byrow = TRUE)
                Measured_D1_Efficacy_SEs <- matrix(Measured_D1_Data$Efficacy_SE, ncol = length(Measured_D1_Data$Efficacy_SE), nrow = n_Simulations, byrow = TRUE)
                SEs_deviated <- (D1_MC_Efficacies[[1]] - Measured_D1_Efficacies) / Measured_D1_Efficacy_SEs
              } else {
                Measured_D1_Data <- Drug1Data[[names(Drug1Data)[i]]][,c("CellLine", "Efficacy", "Efficacy_SE")]
                Measured_D1_Data <- Measured_D1_Data[match(colnames(SEs_deviated), Measured_D1_Data$CellLine),]
                Measured_D1_Efficacies <- matrix(Measured_D1_Data$Efficacy, ncol = length(Measured_D1_Data$Efficacy), nrow = n_Simulations, byrow = TRUE, dimnames = list(NULL, Measured_D1_Data$CellLine))
                Measured_D1_Efficacy_SEs <- matrix(Measured_D1_Data$Efficacy_SE, ncol = length(Measured_D1_Data$Efficacy_SE), nrow = n_Simulations, byrow = TRUE, dimnames = list(NULL, Measured_D1_Data$CellLine))
                D1_MC_Efficacies[[i]] <- Measured_D1_Efficacies + (SEs_deviated * Measured_D1_Efficacy_SEs)
              }
            }
            names(D1_MC_Efficacies) <- names(Drug1Data)
            D2_MC_Efficacies <- as.list(NULL)
            for(i in 1:length(Drug2Data)){
              #Checking if drug 2 is the same as drug 1.
              #If so, not randomly sampling for D2 viabilities. Instead, calculating how many SE's from
              #the measured value each simulated D1 viability fell, and matching that distance in
              #the simulated D2 viabilities. This is done because D1 and D2 viabilities are not independent in
              #this case--they will have come from the same dose-response curve.
                if(Drug1 == Drug2){
                  #Calculating simulated D2 viabilities based on the SE deviations from the simulated D1 viabilities
                    Measured_D2_Data <- Drug2Data[[names(Drug2Data)[i]]][,c("CellLine", "Efficacy", "Efficacy_SE")]
                    Measured_D2_Data <- Measured_D2_Data[match(colnames(SEs_deviated), Measured_D2_Data$CellLine),]
                    Measured_D2_Efficacies <- matrix(Measured_D2_Data$Efficacy, ncol = length(Measured_D2_Data$Efficacy), nrow = n_Simulations, byrow = TRUE, dimnames = list(NULL, Measured_D2_Data$CellLine))
                    Measured_D2_Efficacy_SEs <- matrix(Measured_D2_Data$Efficacy_SE, ncol = length(Measured_D2_Data$Efficacy_SE), nrow = n_Simulations, byrow = TRUE, dimnames = list(NULL, Measured_D2_Data$CellLine))
                    D2_MC_Efficacies[[i]] <- Measured_D2_Efficacies + (SEs_deviated * Measured_D2_Efficacy_SEs)
                } else {
                  #If Drug1 does not equal Drug2, randomly sampling
                    D2_MC_Efficacies[[i]] <- apply(Drug2Data[[i]][,c("Efficacy", "Efficacy_SE")], 1, function(x){rnorm(n = n_Simulations, mean = x[1], sd = x[2])})
                }
            }
            names(D2_MC_Efficacies) <- names(Drug2Data)
            #Cleaning up
              rm(Measured_D1_Data, Measured_D1_Efficacies, Measured_D1_Efficacy_SEs, SEs_deviated)
              if(exists("Measured_D2_Data")){
                rm(Measured_D2_Data, Measured_D2_Efficacies, Measured_D2_Efficacy_SEs)
              }
          #Resampling cell lines for each simulation with replacement to account for the effect of random variation in cell line selection
            Selected_CCLs_Per_Sim <- t(apply(D1_MC_Efficacies[[1]], 1, function(x){return(sample(1:length(x), replace = TRUE))}))
            for(i in 1:length(D1_MC_Efficacies)){
              D1_MC_Efficacies[[i]] <- t(mapply(function(x,y){return(x[y])}, x = split(D1_MC_Efficacies[[i]], row(D1_MC_Efficacies[[i]])), y = split(Selected_CCLs_Per_Sim, row(Selected_CCLs_Per_Sim))))
            }
            for(i in 1:length(D2_MC_Efficacies)){
              D2_MC_Efficacies[[i]] <- t(mapply(function(x,y){return(x[y])}, x = split(D2_MC_Efficacies[[i]], row(D2_MC_Efficacies[[i]])), y = split(Selected_CCLs_Per_Sim, row(Selected_CCLs_Per_Sim))))
            }
          #Looping through all dose comparisons and calculating uncertainties in output values
            for(i in 1:nrow(Dose_Comparisons)){
              #Identifying correct dose data for each drug for this comparison.
                D1_Data <- D1_MC_Efficacies[[which(names(D1_MC_Efficacies) == Dose_Comparisons$Drug1Dose[i])]]
                D2_Data <- D2_MC_Efficacies[[which(names(D2_MC_Efficacies) == Dose_Comparisons$Drug2Dose[i])]]
              #Calculating expected combination efficacy for each cell line
              #using Independent Drug Action
                Combo_Efficacy <- pmax(D1_Data, D2_Data)
              #Calculating average efficacy across sampled cell lines
                D1_efficacies <- rowMeans(D1_Data)
                D2_efficacies <- rowMeans(D2_Data)
                Combo_efficacies <- rowMeans(Combo_Efficacy)
              #Calculating standard errors
                Dose_Comparisons$Mean_Drug1_Efficacy_SE[i] <- sd(D1_efficacies)
                Dose_Comparisons$Mean_Drug2_Efficacy_SE[i] <- sd(D2_efficacies)
                Dose_Comparisons$Mean_Combo_Efficacy_SE[i] <- sd(Combo_efficacies)
              #Calculating 95% confidence intervals for mean treatment efficacies
                index_2.5 <- floor(0.025*n_Simulations)
                index_97.5 <- ceiling(0.975*n_Simulations)
                Dose_Comparisons$`Mean_Drug1_Efficacy_95%_Confidence_Interval`[i] <- paste(sort(D1_efficacies, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                Dose_Comparisons$`Mean_Drug2_Efficacy_95%_Confidence_Interval`[i] <- paste(sort(D2_efficacies, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                Dose_Comparisons$`Mean_Combo_Efficacy_95%_Confidence_Interval`[i] <- paste(sort(Combo_efficacies, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
              #Calculating simulated IDAcomboscores and HRs
                if(Calculate_IDAcomboscore_And_Hazard_Ratio == TRUE){
                  #HRs
                    HRs_vs_D1 <- (1 - Combo_efficacies) / (1 - D1_efficacies)
                    HRs_vs_D2 <- (1 - Combo_efficacies) / (1 - D2_efficacies)
                  #HR standard errors
                    Dose_Comparisons$HR_vs_Drug1_SE[i] <- sd(HRs_vs_D1)
                    Dose_Comparisons$HR_vs_Drug2_SE[i] <- sd(HRs_vs_D2)
                  #HR 95% confidence intervals
                    Dose_Comparisons$`HR_vs_Drug1_95%_Confidence_Interval`[i] <- paste(sort(HRs_vs_D1, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                    Dose_Comparisons$`HR_vs_Drug2_95%_Confidence_Interval`[i] <- paste(sort(HRs_vs_D2, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                  #Calculating HR p values
                    #Drug 1
                      # n_equal_to_1 <- sum(HRs_vs_D1 == 1)
                      # n_less_than_1 <- sum(HRs_vs_D1 < 1)
                      # n_greater_than_1 <- sum(HRs_vs_D1 > 1)
                      # lower_tail <- (n_less_than_1 + 0.5*n_equal_to_1)/length(HRs_vs_D1)
                      # upper_tail <- (n_greater_than_1 + 0.5*n_equal_to_1)/length(HRs_vs_D1)
                      #Calculating one-sided p-value with null hypothesis that HR >= 1
                        Dose_Comparisons$`p_HR_vs_Drug1>=1`[i] <- sum(HRs_vs_D1 >= 1) / length(HRs_vs_D1)
                        #If p-value is 0, setting as p < minimum p value that can be estimated using this many simulations
                          if(Dose_Comparisons$`p_HR_vs_Drug1>=1`[i] == 0){Dose_Comparisons$`p_HR_vs_Drug1>=1`[i] <- paste0("<", 1/n_Simulations)}
                      # #Calculating two-sided p-value with null hypothesis that HR = 1
                      #   Dose_Comparisons$`p_HR_vs_Drug1=1`[i] <- min(lower_tail, upper_tail)*2
                      #   #If p-value is 0, setting as p < minimum p value that can be estimated using this many simulations
                      #     if(Dose_Comparisons$`p_HR_vs_Drug1=1`[i] == 0){Dose_Comparisons$`p_HR_vs_Drug1=1`[i] <- paste0("<", 1/n_Simulations)}
                    #Drug 2
                      # n_equal_to_1 <- sum(HRs_vs_D2 == 1)
                      # n_less_than_1 <- sum(HRs_vs_D2 < 1)
                      # n_greater_than_1 <- sum(HRs_vs_D2 > 1)
                      # lower_tail <- (n_less_than_1 + 0.5*n_equal_to_1)/length(HRs_vs_D2)
                      # upper_tail <- (n_greater_than_1 + 0.5*n_equal_to_1)/length(HRs_vs_D2)
                      #Calculating one-sided p-value with null hypothesis that HR >= 1
                        Dose_Comparisons$`p_HR_vs_Drug2>=1`[i] <- sum(HRs_vs_D2 >= 1) / length(HRs_vs_D2)
                        #If p-value is 0, setting as p < minimum p value that can be estimated using this many simulations
                          if(Dose_Comparisons$`p_HR_vs_Drug2>=1`[i] == 0){Dose_Comparisons$`p_HR_vs_Drug2>=1`[i] <- paste0("<", 1/n_Simulations)}
                      # #Calculating two-sided p-value with null hypothesis that HR = 1
                      #   Dose_Comparisons$`p_HR_vs_Drug2=1`[i] <- min(lower_tail, upper_tail)*2
                      #   #If p-value is 0, setting as p < minimum p value that can be estimated using this many simulations
                      #     if(Dose_Comparisons$`p_HR_vs_Drug2=1`[i] == 0){Dose_Comparisons$`p_HR_vs_Drug2=1`[i] <- paste0("<", 1/n_Simulations)}
                  #IDAcomboscores
                    delta_viabilities <- pmin((1-D1_efficacies), (1-D2_efficacies)) - (1-Combo_efficacies)
                    HR_C_over_Mbests <- pmax(HRs_vs_D1, HRs_vs_D2)
                    MC_IDAcomboscores <- delta_viabilities - delta_viabilities * HR_C_over_Mbests
                  #IDAcomboscore standard error
                    Dose_Comparisons$IDA_Comboscore_SE[i] <- sd(MC_IDAcomboscores)
                  #IDAcomboscore 95% confidence interval
                    Dose_Comparisons$`IDA_Comboscore_95%_Confidence_Interval`[i] <- paste(sort(MC_IDAcomboscores, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                  #Calculating IDAcomboscore p-value
                    Dose_Comparisons$`p_IDA_Comboscore<=0`[i] <- sum(MC_IDAcomboscores <= 0) / length(MC_IDAcomboscores)
                    #If p-value is 0, setting as p < minimum p value that can be estimated using this many simulations
                      if(Dose_Comparisons$`p_IDA_Comboscore<=0`[i] == 0){Dose_Comparisons$`p_IDA_Comboscore<=0`[i] <- paste0("<", 1/n_Simulations)}
                }
            }
          #Cleaning up
            rm(D1_Data, D2_Data, Combo_Efficacy, index_2.5, index_97.5, D1_efficacies, D2_efficacies, Combo_efficacies)
            if(Calculate_IDAcomboscore_And_Hazard_Ratio == TRUE){
              # rm(delta_viabilities, HR_C_over_Mbests, n_equal_to_1, n_less_than_1, n_greater_than_1, lower_tail, upper_tail)
              rm(delta_viabilities, HR_C_over_Mbests, MC_IDAcomboscores, HRs_vs_D1, HRs_vs_D2)
            }
        }
    }

  #Returning Outputs
    #If Calculate_Uncertainty == FALSE, removing SE columns
      if(Calculate_Uncertainty == FALSE){
        Dose_Comparisons <- Dose_Comparisons[,-which(colnames(Dose_Comparisons) %in% c("Mean_Drug1_Efficacy_SE", "Mean_Drug1_Efficacy_95%_Confidence_Interval", "Mean_Drug2_Efficacy_SE", "Mean_Drug2_Efficacy_95%_Confidence_Interval", "Mean_Combo_Efficacy_SE", "Mean_Combo_Efficacy_95%_Confidence_Interval", "HR_vs_Drug1_SE", "HR_vs_Drug1_95%_Confidence_Interval", "p_HR_vs_Drug1>=1", "p_HR_vs_Drug1=1", "HR_vs_Drug2_SE", "HR_vs_Drug2_95%_Confidence_Interval", "p_HR_vs_Drug2>=1", "p_HR_vs_Drug2=1", "IDA_Comboscore_SE", "IDA_Comboscore_95%_Confidence_Interval", "p_IDA_Comboscore<=0"))]
      }
    #If Calculate_IDAcomboscore_And_Hazard_Ratio == FALSE, removing HR and IDAcomboscore columns
      if(Calculate_IDAcomboscore_And_Hazard_Ratio == FALSE){
        Dose_Comparisons <- Dose_Comparisons[,-which(colnames(Dose_Comparisons) %in% c("HR_vs_Drug1", "HR_vs_Drug1_SE", "HR_vs_Drug1_95%_Confidence_Interval", "p_HR_vs_Drug1>=1", "p_HR_vs_Drug1=1", "HR_vs_Drug2", "HR_vs_Drug2_SE", "HR_vs_Drug2_95%_Confidence_Interval", "p_HR_vs_Drug2>=1", "p_HR_vs_Drug2=1", "IDA_Comboscore", "IDA_Comboscore_SE", "IDA_Comboscore_95%_Confidence_Interval", "p_IDA_Comboscore<=0"))]
      }
    #Replacing "Efficacy" with Efficacy_Metric_Name in column names of Dose_Comparisons
      colnames(Dose_Comparisons) <- gsub("Efficacy", Efficacy_Metric_Name, colnames(Dose_Comparisons))
    #Constructing Return_Object
      if(Return_Bootstrap_Values == FALSE){
        Return_Object <- list(Dose_Comparisons, Drug1, Drug2, Usable_CellLines)
        names(Return_Object) <- c("Efficacy_Predictions", "Drug1", "Drug2", "Cell_Lines_Used")
      } else if(Return_Bootstrap_Values == TRUE & Calculate_Uncertainty == TRUE){
        Return_Object <- list(Dose_Comparisons, Drug1, Drug2, Usable_CellLines, D1_MC_Efficacies, D2_MC_Efficacies)
        names(Return_Object) <- c("Efficacy_Predictions", "Drug1", "Drug2", "Cell_Lines_Used", "Bootstrap_Drug1_Efficacies", "Bootstrap_Drug2_Efficacies")
      }
    #Returning output
      return(Return_Object)
}
