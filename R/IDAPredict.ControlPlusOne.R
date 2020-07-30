#' Predicts IDA efficacies for combinations of a control therapy plus a one additional drug
#'
#' This function creates efficacy predictions for combinations of a control therapy + 1 additional drug using monotherapy efficacy data and the assumptions of independent drug action. When data is available for multiple concentrations of the drug to add, efficacy predictions are made for all possible concentration combinations.
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
#' @param Control_Treatment_Drugs A character vector of length > 0 containing the names of the drugs in the control drug treatment for which efficacy predictions are to be made.
#' @param Control_Treatment_Drug_Concentrations A vector of drug concentrations for Control_Treatment_Drugs with the first concentration in Control_Treatment_Drug_Concentrations corresponding to the first drug in Control_Treatment_Drugs etc. Only one concentration may be specified for each drug in the control treatment, but, if a drug is included in both the control and test treatments, there is no need for the same concentration of that drug to be used in both treatments.
#' @param Drug_to_Add A character vector of length 1 containing the name of the drug to add to the control treatment to create a new drug combination for which efficacy predictions are to be made.
#' @param Calculate_Uncertainty A logic vector of length one indicating whether or not a Monte Carlo simulation should be performed to estimate uncertainties in the efficacy predictions based on uncertainties in the monotherapy efficacy measurements. Set TRUE if you wish to calculate uncertainties. Defaults to FALSE.
#' @param Efficacy_SE_Column A character vector of length 1 containing the name of the column in the Monotherapy_Data data frame which contains the standard errors of measured drug efficacies. Must be specified if Calculate_Uncertainty is set to TRUE.
#' @param n_Simulations A positive, integer vector of length 1 with a value >= 40 indicating the number of random samples to be drawn when calculating output efficacy prediction uncertainties. Defaults to 1000.
#' @param Calculate_IDAcomboscore_And_Hazard_Ratio A logic vector of length 1 indicating whether or not IDA-Comboscores and Hazard Ratios (HRs) should be calculated between monotherapies and the drug combination. Set TRUE if so. Should only be set to TRUE for efficacy metrics that range between 0 and 1 (i.e. percent viability). Defaults to FALSE.
#' @param Average_Duplicate_Records A logic vector of length 1 indicating whether or not duplicated records (where a cell line has multiple records for being tested with a given drug at a given concentration) should be averaged. If TRUE, Efficacy values are averaged, and, if Calculate_Uncertainty is also TRUE, Efficacy_SE values are added in quadrature and divided by the number of duplicate records for that cell line/drug/concentration set.
#' @param Return_Bootstrap_Values A logic vector of length 1 indicating whether or not the function should return the Drug1 Efficacies and Drug2 Efficacies simulated in the semi-parametric bootstrap used to estimate the uncertainties of those values. If equal to TRUE, Calculate_Uncertainty must also equal TRUE.
#'
#'@details
#'Uncertainty estimates for values calculated by this function are generated using a semi-parametric bootstrap approach. This is performed in several steps.\enumerate{
#'\item Control treatment efficacies for each drug/concentration are simulated by random sampling from normal distributions with means equal to the provided calculated efficacies and standard deviations equal to the provided efficacy standard errors.
#'\item Drug_to_Add efficacies are simulated in the same fashion as for control treatment efficacies, except in cases when Drug_to_Add is in Control_Treatment_Drugs. In such cases, it is assumed that the efficacy values for Drug_to_Add and its match in Control_Treatment_Drugs are derived from the same dose-response curve, so each simulated efficacy for Drug_to_Add is matched to the corresponding simulated efficacy for that drug in Control_Treatment_Drugs using a standard normal deviate.
#'\item Efficacy predictions are made for the combination of Control_Treatment_Drugs + Drug_to_Add for each cell line and set of simulated efficacies using the assumptions of independent drug action.
#'\item Cell lines are randomly sampled with replacement for each simulation as many times as there are original cell lines. The simulated Control_Treatment_Drugs monotherapy efficacies and combination efficacies are then sampled according to the sampled cell lines for each simulation.
#'\item Mean efficacies are calculated for the control and new combination treatments for each simulation. If specified to do so, these values are then used to calculate simulated HRs and IDAcomboscores.
#'\item The simulated distributions of each efficacy metric are used to estimate uncertainties for those metrics.
#'}
#'
#'@return \itemize{
#'\item If Return_Bootstrap_Values = FALSE, this function returns a list with 4 elements: 1) Either a data frame with the calculated efficacy predictions, or, if an error occurred, a character vector of length one with the error message. 2) A character value with the name of the control therapy 3) A character value with the name of Drug_to_Add 4) A character vector containing the names of the cell lines used to make the efficacy predictions.
#'\item If Return_Bootstrap_Values = TRUE & Calculate_Uncertainty = TRUE, this function returns a list with 6 elements: the first 4 elements are the same as when Return_Bootstrap_Values = FALSE and the fifth and sixth elements are numeric vectors of, respectively, the control therapy and Drug_to_Add viabilities simulated during the semi-parametric bootstrap used to estimate uncertainties.
#'}
#'
#' @examples
#' #Loading Package
#'   library(IDACombo)
#'
#' #Making fake monotherapy dataset
#'   CellLineNames <- rep(c("CL1", "CL2", "CL3", "CL4", "CL5", "CL6"), 6)
#'   DrugNames <- c(rep("D1", 12), rep("D2", 12), rep("D3", 12))
#'   Concentrations <- c(rep(1, 6), rep(2, 6), rep("1.5", 6), rep("3", 6), rep("1.5", 6), rep("3", 6))
#'   Viability <- c(sample(seq(0.4,1,length.out = 10), 6, replace = TRUE),
#'                  sample(seq(0.2,0.8,length.out = 10), 6, replace = TRUE),
#'                  sample(seq(0.4,1,length.out = 10), 6, replace = TRUE),
#'                  sample(seq(0.2,0.6,length.out = 10), 6, replace = TRUE),
#'                  sample(seq(0.38,1,length.out = 10), 6, replace = TRUE),
#'                  sample(seq(0.18,0.6,length.out = 10), 6, replace = TRUE))
#'   Viability_SE <- Viability * sample(seq(0,0.1,length.out = 100), 36, replace = TRUE)
#'   Fake_Data <- data.frame(CellLineNames, DrugNames, Concentrations, Viability, Viability_SE)
#'
#' #Creating efficacy predictions for D1+D2 + D3 without uncertainty calculations
#'   IDAPredict.ControlPlusOne(Monotherapy_Data = Fake_Data,
#'                    Cell_Line_Name_Column = "CellLineNames",
#'                    Drug_Name_Column = "DrugNames",
#'                    Drug_Concentration_Column = "Concentrations",
#'                    Efficacy_Column = "Viability",
#'                    Control_Treatment_Drugs = c("D1", "D2"),
#'                    Control_Treatment_Drug_Concentrations = c(1, 3),
#'                    Drug_to_Add = "D3",
#'                    Calculate_Uncertainty = FALSE,
#'                    LowerEfficacyIsBetterDrugEffect = TRUE,
#'                    Efficacy_Metric_Name = "Viability",
#'                    Calculate_IDAcomboscore_And_Hazard_Ratio = TRUE,
#'                    Average_Duplicate_Records = FALSE)
#'
#' #Creating efficacy predictions for D1+D2 + D3 with uncertainty calculations
#'   IDAPredict.ControlPlusOne(Monotherapy_Data = Fake_Data,
#'                    Cell_Line_Name_Column = "CellLineNames",
#'                    Drug_Name_Column = "DrugNames",
#'                    Drug_Concentration_Column = "Concentrations",
#'                    Efficacy_Column = "Viability",
#'                    Control_Treatment_Drugs = c("D1", "D2"),
#'                    Control_Treatment_Drug_Concentrations = c(1, 3),
#'                    Drug_to_Add = "D3",
#'                    Calculate_Uncertainty = TRUE,
#'                    Efficacy_SE_Column = "Viability_SE",
#'                    LowerEfficacyIsBetterDrugEffect = TRUE,
#'                    Efficacy_Metric_Name = "Viability",
#'                    Calculate_IDAcomboscore_And_Hazard_Ratio = TRUE,
#'                    Average_Duplicate_Records = FALSE)
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
#'   IDAPredict.ControlPlusOne(Monotherapy_Data = Fake_Data,
#'                    Cell_Line_Name_Column = "CellLineNames",
#'                    Drug_Name_Column = "DrugNames",
#'                    Drug_Concentration_Column = "Concentrations",
#'                    Efficacy_Column = "Reduction_in_Viability",
#'                    Control_Treatment_Drugs = c("D1", "D2"),
#'                    Control_Treatment_Drug_Concentrations = c(1, 3),
#'                    Drug_to_Add = "D3",
#'                    Calculate_Uncertainty = TRUE,
#'                    LowerEfficacyIsBetterDrugEffect = FALSE,
#'                    Efficacy_SE_Column = "Reduction_in_Viability_SE",
#'                    Efficacy_Metric_Name = "Reduction_In_Viability",
#'                    Calculate_IDAcomboscore_And_Hazard_Ratio = TRUE,
#'                    Average_Duplicate_Records = FALSE)
#'
#' #Changing efficacy metric to percent growth (range -1 to 1)
#' #Note that calculating Hazard Ratios and IDA-Comboscores is no longer valid, so
#' #Calculate_IDAcomboscore_And_Hazard_Ratio is set to FALSE.
#'   Percent_Growth <- c(sample(seq(0.4,1,length.out = 10), 6, replace = TRUE),
#'                       sample(seq(-0.4,0.2,length.out = 10), 6, replace = TRUE),
#'                       sample(seq(-0.2,0.3,length.out = 10), 6, replace = TRUE),
#'                       sample(seq(-1,0.2,length.out = 10), 6, replace = TRUE),
#'                       sample(seq(-0.22,0.28,length.out = 10), 6, replace = TRUE),
#'                       sample(seq(-1,0.1,length.out = 10), 6, replace = TRUE))
#'   Percent_Growth_SE <- abs(Percent_Growth * sample(seq(0,0.1,length.out = 100), 36, replace = TRUE))
#'   Fake_Data <- data.frame(CellLineNames,
#'                           DrugNames,
#'                           Concentrations,
#'                           Percent_Growth,
#'                           Percent_Growth_SE)
#'   IDAPredict.ControlPlusOne(Monotherapy_Data = Fake_Data,
#'                    Cell_Line_Name_Column = "CellLineNames",
#'                    Drug_Name_Column = "DrugNames",
#'                    Drug_Concentration_Column = "Concentrations",
#'                    Efficacy_Column = "Percent_Growth",
#'                    Control_Treatment_Drugs = c("D1", "D2"),
#'                    Control_Treatment_Drug_Concentrations = c(1, 3),
#'                    Drug_to_Add = "D3",
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
#'                       sample(seq(-1,0.2,length.out = 10), 6, replace = TRUE),
#'                       sample(seq(-0.22,0.28,length.out = 10), 6, replace = TRUE),
#'                       sample(seq(-1,0.1,length.out = 10), 6, replace = TRUE))
#'   Percent_Growth_SE <- abs(Percent_Growth * sample(seq(0,0.1,length.out = 100), 36, replace = TRUE))
#'   Fake_Data_to_add <- data.frame(CellLineNames,
#'                           DrugNames,
#'                           Concentrations,
#'                           Percent_Growth,
#'                           Percent_Growth_SE)
#'   Fake_Data <- rbind(Fake_Data, Fake_Data_to_add)
#'   IDAPredict.ControlPlusOne(Monotherapy_Data = Fake_Data,
#'                    Cell_Line_Name_Column = "CellLineNames",
#'                    Drug_Name_Column = "DrugNames",
#'                    Drug_Concentration_Column = "Concentrations",
#'                    Efficacy_Column = "Percent_Growth",
#'                    Control_Treatment_Drugs = c("D1", "D2"),
#'                    Control_Treatment_Drug_Concentrations = c(1, 3),
#'                    Drug_to_Add = "D3",
#'                    Calculate_Uncertainty = TRUE,
#'                    LowerEfficacyIsBetterDrugEffect = TRUE,
#'                    Efficacy_SE_Column = "Percent_Growth_SE",
#'                    Calculate_IDAcomboscore_And_Hazard_Ratio = FALSE,
#'                    Efficacy_Metric_Name = "Percent_Growth",
#'                    Average_Duplicate_Records = FALSE)
#'
#' #Now setting to average duplicate values.
#'   Fake_Data <- rbind(Fake_Data, Fake_Data_to_add)
#'   IDAPredict.ControlPlusOne(Monotherapy_Data = Fake_Data,
#'                    Cell_Line_Name_Column = "CellLineNames",
#'                    Drug_Name_Column = "DrugNames",
#'                    Drug_Concentration_Column = "Concentrations",
#'                    Efficacy_Column = "Percent_Growth",
#'                    Control_Treatment_Drugs = c("D1", "D2"),
#'                    Control_Treatment_Drug_Concentrations = c(1, 3),
#'                    Drug_to_Add = "D3",
#'                    Calculate_Uncertainty = TRUE,
#'                    LowerEfficacyIsBetterDrugEffect = TRUE,
#'                    Efficacy_SE_Column = "Percent_Growth_SE",
#'                    Calculate_IDAcomboscore_And_Hazard_Ratio = FALSE,
#'                    Efficacy_Metric_Name = "Percent_Growth",
#'                    Average_Duplicate_Records = TRUE)
#' @export
IDAPredict.ControlPlusOne <- function(Monotherapy_Data, Cell_Line_Name_Column, Drug_Name_Column, Drug_Concentration_Column, Efficacy_Column, LowerEfficacyIsBetterDrugEffect, Efficacy_Metric_Name = "Efficacy", Control_Treatment_Drugs, Control_Treatment_Drug_Concentrations, Drug_to_Add, Calculate_Uncertainty = FALSE, Efficacy_SE_Column = NULL, n_Simulations = 1000, Calculate_IDAcomboscore_And_Hazard_Ratio = FALSE, Average_Duplicate_Records = FALSE, Return_Bootstrap_Values = FALSE){

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
    if(! is.vector(Control_Treatment_Drugs) | ! is.character(Control_Treatment_Drugs) | ! length(Control_Treatment_Drugs) > 0){
      stop("Control_Treatment_Drugs is not a character vector with length > 0.")
    }
    if(! all(Control_Treatment_Drugs %in% Monotherapy_Data[,Drug_Name_Column])){
      missing.control.drugs <- Control_Treatment_Drugs[! Control_Treatment_Drugs %in% Monotherapy_Data[,Drug_Name_Column]]
      stop(paste0("No data for the following Control_Treatment_Drugs found in Monotherapy_Data: ", paste(missing.control.drugs, collapse = ", ")))
    }
    if(! is.vector(Control_Treatment_Drug_Concentrations) | ! length(Control_Treatment_Drug_Concentrations) == length(Control_Treatment_Drugs)){
      stop("Control_Treatment_Drug_Concentrations is not a vector with length = length(Control_Treatment_Drugs).")
    }
    if(! is.vector(Drug_to_Add) | ! is.character(Drug_to_Add) | ! length(Drug_to_Add) == 1){
      stop("Drug_to_Add is not a character vector with length 1.")
    }
    if(! Drug_to_Add %in% Monotherapy_Data[,Drug_Name_Column]){
      stop(paste0("No ", Drug_to_Add, " data found in Monotherapy_Data."))
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
      if(! n_Simulations%%1==0 | ! n_Simulations > 0){
        stop("Calculate_Uncertainty is TRUE, but n_Simulations is not a positive, non-zero integer.")
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


  #Creating single name for control treatment
    Control_Treatment_Name <- paste0("combo:", paste(Control_Treatment_Drugs, collapse = "+"))

  #Organizing data into standard format based on column names provided for each desired set of information
  #Also subsetting to only include data pertaining to Control_Treatment_Drugs and Drug_to_Add
    if(Calculate_Uncertainty == TRUE){
      Data <- Monotherapy_Data[Monotherapy_Data[,Drug_Name_Column] %in% c(Control_Treatment_Drugs, Drug_to_Add),c(Cell_Line_Name_Column, Drug_Name_Column, Drug_Concentration_Column, Efficacy_Column, Efficacy_SE_Column)]
      colnames(Data) <- c("CellLine", "Drug", "Conc", "Efficacy", "Efficacy_SE")
    } else {
      Data <- Monotherapy_Data[Monotherapy_Data[,Drug_Name_Column] %in% c(Control_Treatment_Drugs, Drug_to_Add),c(Cell_Line_Name_Column, Drug_Name_Column, Drug_Concentration_Column, Efficacy_Column)]
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

  #Subsetting into control group data with specified control concentrations
    ControlData <- list(NULL)
    for(i in 1:length(Control_Treatment_Drugs)){
      ControlData[[i]] <- Data[Data$Drug %in% Control_Treatment_Drugs[i],]
      #Checking that all provided concentrations are available for their respective drugs
        if(! Control_Treatment_Drug_Concentrations[i] %in% ControlData[[i]]$Conc){
          stop(paste0(Control_Treatment_Drug_Concentrations[i], " concentration is unavailable for ", Control_Treatment_Drugs[i], " in control treatment."))
        } else {
          ControlData[[i]] <- ControlData[[i]][ControlData[[i]]$Conc %in% Control_Treatment_Drug_Concentrations[i],]
        }
    }
    names(ControlData) <- Control_Treatment_Drugs

  #Collapsing control treatment into "single drug" treatment which represents activity of the control combination
    #Finding cell line overlap between all drugs in control treatment
      ControlCellLines <- list(NULL)
      for(i in 1:length(ControlData)){
        ControlCellLines[[i]] <- sort(unique(ControlData[[i]]$CellLine))
      }
      Usable_CellLines <- sort(unique(unlist(ControlCellLines)))
      for(i in 1:length(ControlCellLines)){
        Usable_CellLines <- Usable_CellLines[Usable_CellLines %in% ControlCellLines[[i]]]
      }
      rm(ControlCellLines)
    #Checking if at least 2 cell lines remain for all control drugs. If not, exiting with no
    #predictions and a warning.
      if(! length(Usable_CellLines) >= 2){
        #Returning NA predictions with warning due to too few cell lines.
        warning(paste0("<2 overlapping cell lines available for combination of ", gsub("^combo:", "", Control_Treatment_Name), " + ", Drug_to_Add, "."))
        if(Return_Bootstrap_Values == FALSE){
          Return_Object <- list("Less than 2 overlapping cell lines available.", gsub("^combo:", "", Control_Treatment_Name), Drug_to_Add, Usable_CellLines)
          names(Return_Object) <- c("Efficacy_Predictions", "Control_Treatment", "Drug_to_Add", "Cell_Lines_Used")
          return(Return_Object)
        } else if(Return_Bootstrap_Values == TRUE){
          Return_Object <- list("Less than 2 overlapping cell lines available.", gsub("^combo:", "", Control_Treatment_Name), Drug_to_Add, Usable_CellLines, NULL, NULL)
          names(Return_Object) <- c("Efficacy_Predictions", "Control_Treatment", "Drug_to_Add", "Cell_Lines_Used", "Bootstrap_Control_Efficacies", "Boostrap_Drug_to_Add_Efficacies")
          return(Return_Object)
        }
      }

    #Subsetting drug data to only include overlapping cell lines
      for(i in 1:length(ControlData)){
        ControlData[[i]] <- ControlData[[i]][ControlData[[i]]$CellLine %in% Usable_CellLines,]
      }
    #Checking that cell lines aren't duplicated in each control drug dataset
    #If Average_Duplicate_Records == FALSE, removing cell line duplicates with warning if duplicates are found.
    #If Average_Duplicate_Records == TRUE, averaging duplicate records without warning.
      for(i in 1:length(Control_Treatment_Drugs)){
        CL_Conc <- paste(ControlData[[i]]$CellLine, ControlData[[i]]$Conc, sep = "_")
        Dups <- CL_Conc[duplicated(CL_Conc)]
        if(length(Dups) > 0 & Average_Duplicate_Records == FALSE){
          warning(paste0("Duplicated information found for the following cell lines and ", Control_Treatment_Drugs[i], " concentrations in the control treatment. Average_Duplicate_Records = FALSE so duplicates removed: ", paste(Dups, collapse = ", ")))
          ControlData[[i]] <- ControlData[[i]][! duplicated(CL_Conc),]
        } else if(length(Dups) > 0 & Average_Duplicate_Records == TRUE){
          if(Calculate_Uncertainty == FALSE){
            #Simply averaging efficacy values
              colnames <- colnames(ControlData[[i]])
              ControlData[[i]] <- aggregate(ControlData[[i]]$Efficacy, by = list(ControlData[[i]]$CellLine, ControlData[[i]]$Drug, ControlData[[i]]$Conc), FUN = mean)
              colnames(ControlData[[i]]) <- colnames
          } else if(Calculate_Uncertainty == TRUE){
            #Averaging efficacy
              Efficacy_Average <- aggregate(ControlData[[i]]$Efficacy, by = list(ControlData[[i]]$CellLine, ControlData[[i]]$Drug, ControlData[[i]]$Conc), FUN = mean)
              colnames(Efficacy_Average) <- c("CellLine", "Drug", "Conc", "Efficacy")
            #Calculating uncertainty in averaged efficacy by adding efficacy uncertainties in quadrature and dividing by number of values used in average
              Efficacy_Average_Uncertainties <- aggregate(ControlData[[i]]$Efficacy_SE, by = list(ControlData[[i]]$CellLine, ControlData[[i]]$Drug, ControlData[[i]]$Conc), FUN = function(x){sqrt(sum(x^2))/length(x)})
              colnames(Efficacy_Average_Uncertainties) <- c("CellLine", "Drug", "Conc", "Efficacy_SE")
            #Combining results
              ControlData[[i]] <- merge(Efficacy_Average, Efficacy_Average_Uncertainties)
          }
        }
      }
    #Ordering cell lines the same for each drug
      for(i in 1:length(ControlData)){
        ControlData[[i]] <- ControlData[[i]][order(ControlData[[i]]$Conc, ControlData[[i]]$CellLine),]
      }

  #Subsetting Drug_to_Add data
    Drug_to_AddData <- Data[Data$Drug %in% Drug_to_Add,]
    rm(Data)

  #Finding cell line overlap between all drugs
    Control_Cell_Lines <- sort(unique(ControlData[[1]]$CellLine))
    Usable_CellLines <- sort(unique(Control_Cell_Lines[Control_Cell_Lines %in% Drug_to_AddData$CellLine]))

  #Checking again if at least 2 cell lines remain for all control drugs. If not, exiting with no
  #predictions and a warning.
    if(! length(Usable_CellLines) >= 2){
      #Returning NA predictions with warning due to too few cell lines.
      warning(paste0("<2 overlapping cell lines available for combination of ", gsub("^combo:", "", Control_Treatment_Name), " + ", Drug_to_Add, "."))
      if(Return_Bootstrap_Values == FALSE){
        Return_Object <- list("Less than 2 overlapping cell lines available.", gsub("^combo:", "", Control_Treatment_Name), Drug_to_Add, Usable_CellLines)
        names(Return_Object) <- c("Efficacy_Predictions", "Control_Treatment", "Drug_to_Add", "Cell_Lines_Used")
        return(Return_Object)
      } else if(Return_Bootstrap_Values == TRUE){
        Return_Object <- list("Less than 2 overlapping cell lines available.", gsub("^combo:", "", Control_Treatment_Name), Drug_to_Add, Usable_CellLines, NULL, NULL)
        names(Return_Object) <- c("Efficacy_Predictions", "Control_Treatment", "Drug_to_Add", "Cell_Lines_Used", "Bootstrap_Control_Efficacies", "Boostrap_Drug_to_Add_Efficacies")
        return(Return_Object)
      }
    }

  #Subsetting drug data to only include overlapping cell lines and ordering same
  #for each drug.
    for(i in 1:length(ControlData)){
      ControlData[[i]] <- ControlData[[i]][ControlData[[i]]$CellLine %in% Usable_CellLines,]
      ControlData[[i]] <- ControlData[[i]][order(ControlData[[i]]$Conc, ControlData[[i]]$CellLine),]
    }
    Drug_to_AddData <- Drug_to_AddData[Drug_to_AddData$CellLine %in% Usable_CellLines,]
    Drug_to_AddData <- Drug_to_AddData[order(Drug_to_AddData$Conc, Drug_to_AddData$CellLine),]

  #Checking that cell lines aren't duplicated in Drug_to_AddData for a given drug dose.
  #If Average_Duplicate_Records == FALSE, removing cell line duplicates with warning if duplicates are found.
  #If Average_Duplicate_Records == TRUE, averaging duplicate records without warning.
    D2_CL_Conc <- paste(Drug_to_AddData$CellLine, Drug_to_AddData$Conc, sep = "_")
    D2Dups <- D2_CL_Conc[duplicated(D2_CL_Conc)]
    if(length(D2Dups) > 0 & Average_Duplicate_Records == FALSE){
      warning(paste0("Duplicated information found for the following cell lines and ", Drug_to_Add, " concentrations. Average_Duplicate_Records = FALSE so duplicates removed: ", paste(D2Dups, collapse = ", ")))
      Drug_to_AddData <- Drug_to_AddData[! duplicated(D2_CL_Conc),]
    } else if(length(D2Dups) > 0 & Average_Duplicate_Records == TRUE){
      if(Calculate_Uncertainty == FALSE){
        #Simply averaging efficacy values
          d2.colnames <- colnames(Drug_to_AddData)
          Drug_to_AddData <- aggregate(Drug_to_AddData$Efficacy, by = list(Drug_to_AddData$CellLine, Drug_to_AddData$Drug, Drug_to_AddData$Conc), FUN = mean)
          colnames(Drug_to_AddData) <- d2.colnames
      } else if(Calculate_Uncertainty == TRUE){
        #Averaging efficacy
          d2.Efficacy_Average <- aggregate(Drug_to_AddData$Efficacy, by = list(Drug_to_AddData$CellLine, Drug_to_AddData$Drug, Drug_to_AddData$Conc), FUN = mean)
          colnames(d2.Efficacy_Average) <- c("CellLine", "Drug", "Conc", "Efficacy")
        #Calculating uncertainty in averaged efficacy by adding efficacy uncertainties in quadrature and dividing by number of values used in average
          d2.Efficacy_Average_Uncertainties <- aggregate(Drug_to_AddData$Efficacy_SE, by = list(Drug_to_AddData$CellLine, Drug_to_AddData$Drug, Drug_to_AddData$Conc), FUN = function(x){sqrt(sum(x^2))/length(x)})
          colnames(d2.Efficacy_Average_Uncertainties) <- c("CellLine", "Drug", "Conc", "Efficacy_SE")
        #Combining results
          Drug_to_AddData <- merge(d2.Efficacy_Average, d2.Efficacy_Average_Uncertainties)
      }
    }

  #Checking that all drug concentrations for each drug are available for each cell line
  #Omitting cell lines that are missing concentrations for 1 or more drugs
    AllData <- append(ControlData, list(Drug_to_AddData))
    names(AllData)[length(AllData)] <- Drug_to_Add
    CLs_per_dose <- sapply(AllData, function(x){as.data.frame.table(table(x$Conc), stringsAsFactors = FALSE)[,2]})
    if(! all(unlist(unique(c(CLs_per_dose))) == length(Usable_CellLines))){
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
        warning(paste0("The following cell lines are missing efficacy information for one or more drug concentrations in the combination of ", gsub("^combo:", "", Control_Treatment_Name), " + ", Drug_to_Add, ". These cell lines have been omitted from the analysis: ", paste(CLs_missing_doses, collapse = ", ")))
      #Re-subsetting drug data to only include overlapping cell lines and ordering same
      #for each drug with cell lines that had missing information removed.
        Usable_CellLines <- Usable_CellLines[! Usable_CellLines %in% CLs_missing_doses]
        for(i in 1:length(ControlData)){
          ControlData[[i]] <- ControlData[[i]][ControlData[[i]]$CellLine %in% Usable_CellLines,]
          ControlData[[i]] <- ControlData[[i]][order(ControlData[[i]]$Conc, ControlData[[i]]$CellLine),]
        }
        Drug_to_AddData <- Drug_to_AddData[Drug_to_AddData$CellLine %in% Usable_CellLines,]
        Drug_to_AddData <- Drug_to_AddData[order(Drug_to_AddData$Conc, Drug_to_AddData$CellLine),]
        rm(AllData, CLs_per_dose, n_doses, CLs_missing_doses, CL_dose_count)
    } else {
      rm(AllData, CLs_per_dose)
    }

  #Checking again if at least 2 cell lines remain for all control drugs. If not, exiting with no
  #predictions and a warning.
    if(! length(Usable_CellLines) >= 2){
      #Returning NA predictions with warning due to too few cell lines.
      warning(paste0("<2 overlapping cell lines available for combination of ", gsub("^combo:", "", Control_Treatment_Name), " + ", Drug_to_Add, "."))
      if(Return_Bootstrap_Values == FALSE){
        Return_Object <- list("Less than 2 overlapping cell lines available.", gsub("^combo:", "", Control_Treatment_Name), Drug_to_Add, Usable_CellLines)
        names(Return_Object) <- c("Efficacy_Predictions", "Control_Treatment", "Drug_to_Add", "Cell_Lines_Used")
        return(Return_Object)
      } else if(Return_Bootstrap_Values == TRUE){
        Return_Object <- list("Less than 2 overlapping cell lines available.", gsub("^combo:", "", Control_Treatment_Name), Drug_to_Add, Usable_CellLines, NULL, NULL)
        names(Return_Object) <- c("Efficacy_Predictions", "Control_Treatment", "Drug_to_Add", "Cell_Lines_Used", "Bootstrap_Control_Efficacies", "Boostrap_Drug_to_Add_Efficacies")
        return(Return_Object)
      }
    }

  #Dividing Drug_to_AddData by drug concentration
    Drug_to_AddData <- split(Drug_to_AddData, as.factor(Drug_to_AddData$Conc))

  #Identifying all concentration comparisons that need to be made.
    D1_concentrations <- paste(Control_Treatment_Drug_Concentrations, collapse = ",")
    D2_concentrations <- names(Drug_to_AddData)
    Dose_Comparisons <- expand.grid(D1_concentrations, D2_concentrations, stringsAsFactors = FALSE)
    colnames(Dose_Comparisons) <- c("Control_Treatment_Doses", "Drug_to_Add_Dose")

  #Adding extra columns to Dose_Comparisons to store predicted efficacy results
    Dose_Comparisons$Mean_Control_Treatment_Efficacy <- NA
    Dose_Comparisons$Mean_Control_Treatment_Efficacy_SE <- NA
    Dose_Comparisons$`Mean_Control_Treatment_Efficacy_95%_Confidence_Interval` <- NA
    Dose_Comparisons$Mean_Drug_to_Add_Efficacy <- NA
    Dose_Comparisons$Mean_Drug_to_Add_Efficacy_SE <- NA
    Dose_Comparisons$`Mean_Drug_to_Add_Efficacy_95%_Confidence_Interval` <- NA
    Dose_Comparisons$Mean_Combo_Efficacy <- NA
    Dose_Comparisons$Mean_Combo_Efficacy_SE <- NA
    Dose_Comparisons$`Mean_Combo_Efficacy_95%_Confidence_Interval` <- NA
    Dose_Comparisons$HR_vs_Control_Treatment <- NA
    Dose_Comparisons$HR_vs_Control_Treatment_SE <- NA
    Dose_Comparisons$`HR_vs_Control_Treatment_95%_Confidence_Interval` <- NA
    Dose_Comparisons$`p_HR_vs_Control_Treatment>=1` <- NA
    Dose_Comparisons$HR_vs_Drug_to_Add <- NA
    Dose_Comparisons$HR_vs_Drug_to_Add_SE <- NA
    Dose_Comparisons$`HR_vs_Drug_to_Add_95%_Confidence_Interval` <- NA
    Dose_Comparisons$`p_HR_vs_Drug_to_Add>=1` <- NA
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
            D1_Data_Efficacy <- do.call(pmin, lapply(ControlData, function(x){return(x$Efficacy)}))
            D2_Data <- Drug_to_AddData[[which(names(Drug_to_AddData) == Dose_Comparisons$Drug_to_Add_Dose[i])]]
          #Calculating expected combination efficacy for each cell line
          #using Independent Drug Action
            Combo_Efficacy <- pmin(D1_Data_Efficacy, D2_Data$Efficacy)
          #Calculating average efficacy across all cell lines
            Dose_Comparisons$Mean_Control_Treatment_Efficacy[i] <- mean(D1_Data_Efficacy)
            Dose_Comparisons$Mean_Drug_to_Add_Efficacy[i] <- mean(D2_Data$Efficacy)
            Dose_Comparisons$Mean_Combo_Efficacy[i] <- mean(Combo_Efficacy)
            if(Calculate_IDAcomboscore_And_Hazard_Ratio == TRUE){
              Dose_Comparisons$HR_vs_Control_Treatment[i] <- Dose_Comparisons$Mean_Combo_Efficacy[i] / Dose_Comparisons$Mean_Control_Treatment_Efficacy[i]
              Dose_Comparisons$HR_vs_Drug_to_Add[i] <- Dose_Comparisons$Mean_Combo_Efficacy[i] / Dose_Comparisons$Mean_Drug_to_Add_Efficacy[i]
              delta_viability <- min(Dose_Comparisons$Mean_Control_Treatment_Efficacy[i], Dose_Comparisons$Mean_Drug_to_Add_Efficacy[i]) - Dose_Comparisons$Mean_Combo_Efficacy[i]
              HR_C_over_Mbest <- max(Dose_Comparisons$HR_vs_Control_Treatment[i], Dose_Comparisons$HR_vs_Drug_to_Add[i])
              Dose_Comparisons$IDA_Comboscore[i] <- delta_viability - delta_viability * HR_C_over_Mbest
            }
          }
        rm(D1_Data_Efficacy, D2_Data, Combo_Efficacy)
        if(Calculate_IDAcomboscore_And_Hazard_Ratio == TRUE){
          rm(delta_viability, HR_C_over_Mbest)
        }

      #If Calculate_Uncertainty == TRUE, doing Monte Carlo simulation to estimate
      #uncertainties in output parameters based on uncertainties in monotherapy efficacies.
        if(Calculate_Uncertainty == TRUE){
          #Looping through each drug for each treatment and simulating efficacies based on
          #measured efficacies and SE's
            D1_MC_Efficacies <- as.list(NULL)
            SEs_deviated_list <- as.list(NULL)
            for(i in 1:length(ControlData)){
                D1_MC_Efficacies[[i]] <- apply(ControlData[[i]][,c("Efficacy", "Efficacy_SE")], 1, function(x){rnorm(n = n_Simulations, mean = x[1], sd = x[2])})
                colnames(D1_MC_Efficacies[[i]]) <- ControlData[[i]]$CellLine
                Measured_D1_Data <- ControlData[[names(ControlData)[i]]][,c("CellLine", "Efficacy", "Efficacy_SE")]
                Measured_D1_Data <- Measured_D1_Data[match(colnames(D1_MC_Efficacies[[i]]), Measured_D1_Data$CellLine),]
                Measured_D1_Efficacies <- matrix(Measured_D1_Data$Efficacy, ncol = length(Measured_D1_Data$Efficacy), nrow = n_Simulations, byrow = TRUE)
                Measured_D1_Efficacy_SEs <- matrix(Measured_D1_Data$Efficacy_SE, ncol = length(Measured_D1_Data$Efficacy_SE), nrow = n_Simulations, byrow = TRUE)
                SEs_deviated_list[[i]] <- (D1_MC_Efficacies[[i]] - Measured_D1_Efficacies) / Measured_D1_Efficacy_SEs
            }
            names(SEs_deviated_list) <- names(ControlData)
            names(D1_MC_Efficacies) <- names(ControlData)
            D2_MC_Efficacies <- as.list(NULL)
            for(i in 1:length(Drug_to_AddData)){
              #Checking if Drug_to_Add is one of the control treatment drugs.
              #If so, not randomly sampling for Drug_to_Add viabilities. Instead, calculating how many SE's from
              #the measured value each simulated control drug viability fell, and matching that distance in
              #the simulated Drug_to_Add viabilities. This is done because Drug_to_Add viabilities are not independent in
              #this case--they will have come from the same dose-response curve as was used for the matching control treatment drug.
                if(Drug_to_Add %in% names(ControlData)){
                  #Finding the number of SEs away from the measured value each simulated value is for the matching control drug
                    SEs_deviated <- SEs_deviated_list[[Drug_to_Add]]
                  #Calculating simulated D2 viabilities based on the SE deviations from the simulated D1 viabilities
                    Measured_D2_Data <- Drug_to_AddData[[names(Drug_to_AddData)[i]]][,c("CellLine", "Efficacy", "Efficacy_SE")]
                    Measured_D2_Data <- Measured_D2_Data[match(colnames(SEs_deviated), Measured_D2_Data$CellLine),]
                    Measured_D2_Efficacies <- matrix(Measured_D2_Data$Efficacy, ncol = length(Measured_D2_Data$Efficacy), nrow = n_Simulations, byrow = TRUE, dimnames = list(NULL, Measured_D2_Data$CellLine))
                    Measured_D2_Efficacy_SEs <- matrix(Measured_D2_Data$Efficacy_SE, ncol = length(Measured_D2_Data$Efficacy_SE), nrow = n_Simulations, byrow = TRUE, dimnames = list(NULL, Measured_D2_Data$CellLine))
                    D2_MC_Efficacies[[i]] <- Measured_D2_Efficacies + (SEs_deviated * Measured_D2_Efficacy_SEs)
                } else {
                  #If D2 drug is not in D1 therapy, randomly sampling
                    D2_MC_Efficacies[[i]] <- apply(Drug_to_AddData[[i]][,c("Efficacy", "Efficacy_SE")], 1, function(x){rnorm(n = n_Simulations, mean = x[1], sd = x[2])})
                }
            }
            names(D2_MC_Efficacies) <- names(Drug_to_AddData)
            #Cleaning up
              rm(Measured_D1_Data, Measured_D1_Efficacies, Measured_D1_Efficacy_SEs, SEs_deviated_list)
              if(exists("Measured_D2_Data")){
                rm(SEs_deviated, Measured_D2_Data, Measured_D2_Efficacies, Measured_D2_Efficacy_SEs)
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
                D1_Data <- do.call(pmin, D1_MC_Efficacies)
                D2_Data <- D2_MC_Efficacies[[which(names(D2_MC_Efficacies) == Dose_Comparisons$Drug_to_Add_Dose[i])]]
              #Calculating expected combination efficacy for each cell line
              #using Independent Drug Action
                Combo_Efficacy <- pmin(D1_Data, D2_Data)
              #Calculating average efficacy across sampled cell lines
                D1_efficacies <- rowMeans(D1_Data)
                D2_efficacies <- rowMeans(D2_Data)
                Combo_efficacies <- rowMeans(Combo_Efficacy)
              #Calculating standard errors
                Dose_Comparisons$Mean_Control_Treatment_Efficacy_SE[i] <- sd(D1_efficacies)
                Dose_Comparisons$Mean_Drug_to_Add_Efficacy_SE[i] <- sd(D2_efficacies)
                Dose_Comparisons$Mean_Combo_Efficacy_SE[i] <- sd(Combo_efficacies)
              #Calculating 95% confidence intervals for mean treatment efficacies
                index_2.5 <- floor(0.025*n_Simulations)
                index_97.5 <- ceiling(0.975*n_Simulations)
                Dose_Comparisons$`Mean_Control_Treatment_Efficacy_95%_Confidence_Interval`[i] <- paste(sort(D1_efficacies, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                Dose_Comparisons$`Mean_Drug_to_Add_Efficacy_95%_Confidence_Interval`[i] <- paste(sort(D2_efficacies, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                Dose_Comparisons$`Mean_Combo_Efficacy_95%_Confidence_Interval`[i] <- paste(sort(Combo_efficacies, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
              #Calculating simulated IDAcomboscores and HRs
                if(Calculate_IDAcomboscore_And_Hazard_Ratio == TRUE){
                  #HRs
                    HRs_vs_D1 <- Combo_efficacies / D1_efficacies
                    HRs_vs_D2 <- Combo_efficacies / D2_efficacies
                  #HR standard errors
                    Dose_Comparisons$HR_vs_Control_Treatment_SE[i] <- sd(HRs_vs_D1)
                    Dose_Comparisons$HR_vs_Drug_to_Add_SE[i] <- sd(HRs_vs_D2)
                  #HR 95% confidence intervals
                    Dose_Comparisons$`HR_vs_Control_Treatment_95%_Confidence_Interval`[i] <- paste(sort(HRs_vs_D1, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                    Dose_Comparisons$`HR_vs_Drug_to_Add_95%_Confidence_Interval`[i] <- paste(sort(HRs_vs_D2, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                  #Calculating HR p values
                    #Control therapy
                      #Calculating one-sided p-value with null hypothesis that HR >= 1
                        Dose_Comparisons$`p_HR_vs_Control_Treatment>=1`[i] <- sum(HRs_vs_D1 >= 1) / length(HRs_vs_D1)
                        #If p-value is 0, setting as p < minimum p value that can be estimated using this many simulations
                          if(Dose_Comparisons$`p_HR_vs_Control_Treatment>=1`[i] == 0){Dose_Comparisons$`p_HR_vs_Control_Treatment>=1`[i] <- paste0("<", 1/n_Simulations)}
                    #Drug_to_Add
                      #Calculating one-sided p-value with null hypothesis that HR >= 1
                        Dose_Comparisons$`p_HR_vs_Drug_to_Add>=1`[i] <- sum(HRs_vs_D2 >= 1) / length(HRs_vs_D2)
                        #If p-value is 0, setting as p < minimum p value that can be estimated using this many simulations
                          if(Dose_Comparisons$`p_HR_vs_Drug_to_Add>=1`[i] == 0){Dose_Comparisons$`p_HR_vs_Drug_to_Add>=1`[i] <- paste0("<", 1/n_Simulations)}
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
            D1_Data_Efficacy <- do.call(pmax, lapply(ControlData, function(x){return(x$Efficacy)}))
            D2_Data <- Drug_to_AddData[[which(names(Drug_to_AddData) == Dose_Comparisons$Drug_to_Add_Dose[i])]]
          #Calculating expected combination efficacy for each cell line
          #using Independent Drug Action
            Combo_Efficacy <- pmax(D1_Data_Efficacy, D2_Data$Efficacy)
          #Calculating average efficacy across all cell lines
            Dose_Comparisons$Mean_Control_Treatment_Efficacy[i] <- mean(D1_Data_Efficacy)
            Dose_Comparisons$Mean_Drug_to_Add_Efficacy[i] <- mean(D2_Data$Efficacy)
            Dose_Comparisons$Mean_Combo_Efficacy[i] <- mean(Combo_Efficacy)
            if(Calculate_IDAcomboscore_And_Hazard_Ratio == TRUE){
              Dose_Comparisons$HR_vs_Control_Treatment[i] <- (1 - Dose_Comparisons$Mean_Combo_Efficacy[i]) / (1 - Dose_Comparisons$Mean_Control_Treatment_Efficacy[i])
              Dose_Comparisons$HR_vs_Drug_to_Add[i] <- (1 - Dose_Comparisons$Mean_Combo_Efficacy[i]) / (1 - Dose_Comparisons$Mean_Drug_to_Add_Efficacy[i])
              delta_viability <- min((1 - Dose_Comparisons$Mean_Control_Treatment_Efficacy[i]), (1 - Dose_Comparisons$Mean_Drug_to_Add_Efficacy[i])) - (1 - Dose_Comparisons$Mean_Combo_Efficacy[i])
              HR_C_over_Mbest <- max(Dose_Comparisons$HR_vs_Control_Treatment[i], Dose_Comparisons$HR_vs_Drug_to_Add[i])
              Dose_Comparisons$IDA_Comboscore[i] <- delta_viability - delta_viability * HR_C_over_Mbest
            }
          }
        rm(D1_Data_Efficacy, D2_Data, Combo_Efficacy)
        if(Calculate_IDAcomboscore_And_Hazard_Ratio == TRUE){
          rm(delta_viability, HR_C_over_Mbest)
        }

      #If Calculate_Uncertainty == TRUE, doing Monte Carlo simulation to estimate
      #uncertainties in output parameters based on uncertainties in monotherapy efficacies.
        if(Calculate_Uncertainty == TRUE){
          #Looping through each drug for each treatment and simulating efficacies based on
          #measured efficacies and SE's
            D1_MC_Efficacies <- as.list(NULL)
            SEs_deviated_list <- as.list(NULL)
            for(i in 1:length(ControlData)){
                D1_MC_Efficacies[[i]] <- apply(ControlData[[i]][,c("Efficacy", "Efficacy_SE")], 1, function(x){rnorm(n = n_Simulations, mean = x[1], sd = x[2])})
                colnames(D1_MC_Efficacies[[i]]) <- ControlData[[i]]$CellLine
                Measured_D1_Data <- ControlData[[names(ControlData)[i]]][,c("CellLine", "Efficacy", "Efficacy_SE")]
                Measured_D1_Data <- Measured_D1_Data[match(colnames(D1_MC_Efficacies[[i]]), Measured_D1_Data$CellLine),]
                Measured_D1_Efficacies <- matrix(Measured_D1_Data$Efficacy, ncol = length(Measured_D1_Data$Efficacy), nrow = n_Simulations, byrow = TRUE)
                Measured_D1_Efficacy_SEs <- matrix(Measured_D1_Data$Efficacy_SE, ncol = length(Measured_D1_Data$Efficacy_SE), nrow = n_Simulations, byrow = TRUE)
                SEs_deviated_list[[i]] <- (D1_MC_Efficacies[[i]] - Measured_D1_Efficacies) / Measured_D1_Efficacy_SEs
            }
            names(SEs_deviated_list) <- names(ControlData)
            names(D1_MC_Efficacies) <- names(ControlData)
            D2_MC_Efficacies <- as.list(NULL)
            for(i in 1:length(Drug_to_AddData)){
              #Checking if Drug_to_Add is one of the control treatment drugs.
              #If so, not randomly sampling for Drug_to_Add viabilities. Instead, calculating how many SE's from
              #the measured value each simulated control drug viability fell, and matching that distance in
              #the simulated Drug_to_Add viabilities. This is done because Drug_to_Add viabilities are not independent in
              #this case--they will have come from the same dose-response curve as was used for the matching control treatment drug.
                if(Drug_to_Add %in% names(ControlData)){
                  #Finding the number of SEs away from the measured value each simulated value is for the matching control drug
                    SEs_deviated <- SEs_deviated_list[[Drug_to_Add]]
                  #Calculating simulated D2 viabilities based on the SE deviations from the simulated D1 viabilities
                    Measured_D2_Data <- Drug_to_AddData[[names(Drug_to_AddData)[i]]][,c("CellLine", "Efficacy", "Efficacy_SE")]
                    Measured_D2_Data <- Measured_D2_Data[match(colnames(SEs_deviated), Measured_D2_Data$CellLine),]
                    Measured_D2_Efficacies <- matrix(Measured_D2_Data$Efficacy, ncol = length(Measured_D2_Data$Efficacy), nrow = n_Simulations, byrow = TRUE, dimnames = list(NULL, Measured_D2_Data$CellLine))
                    Measured_D2_Efficacy_SEs <- matrix(Measured_D2_Data$Efficacy_SE, ncol = length(Measured_D2_Data$Efficacy_SE), nrow = n_Simulations, byrow = TRUE, dimnames = list(NULL, Measured_D2_Data$CellLine))
                    D2_MC_Efficacies[[i]] <- Measured_D2_Efficacies + (SEs_deviated * Measured_D2_Efficacy_SEs)
                } else {
                  #If D2 drug is not in D1 therapy, randomly sampling
                    D2_MC_Efficacies[[i]] <- apply(Drug_to_AddData[[i]][,c("Efficacy", "Efficacy_SE")], 1, function(x){rnorm(n = n_Simulations, mean = x[1], sd = x[2])})
                }
            }
            names(D2_MC_Efficacies) <- names(Drug_to_AddData)
            #Cleaning up
              rm(Measured_D1_Data, Measured_D1_Efficacies, Measured_D1_Efficacy_SEs, SEs_deviated_list)
              if(exists("Measured_D2_Data")){
                rm(SEs_deviated, Measured_D2_Data, Measured_D2_Efficacies, Measured_D2_Efficacy_SEs)
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
                D1_Data <- do.call(pmax, D1_MC_Efficacies)
                D2_Data <- D2_MC_Efficacies[[which(names(D2_MC_Efficacies) == Dose_Comparisons$Drug_to_Add_Dose[i])]]
              #Calculating expected combination efficacy for each cell line
              #using Independent Drug Action
                Combo_Efficacy <- pmax(D1_Data, D2_Data)
              #Calculating average efficacy across sampled cell lines
                D1_efficacies <- rowMeans(D1_Data)
                D2_efficacies <- rowMeans(D2_Data)
                Combo_efficacies <- rowMeans(Combo_Efficacy)
              #Calculating standard errors
                Dose_Comparisons$Mean_Control_Treatment_Efficacy_SE[i] <- sd(D1_efficacies)
                Dose_Comparisons$Mean_Drug_to_Add_Efficacy_SE[i] <- sd(D2_efficacies)
                Dose_Comparisons$Mean_Combo_Efficacy_SE[i] <- sd(Combo_efficacies)
              #Calculating 95% confidence intervals for mean treatment efficacies
                index_2.5 <- floor(0.025*n_Simulations)
                index_97.5 <- ceiling(0.975*n_Simulations)
                Dose_Comparisons$`Mean_Control_Treatment_Efficacy_95%_Confidence_Interval`[i] <- paste(sort(D1_efficacies, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                Dose_Comparisons$`Mean_Drug_to_Add_Efficacy_95%_Confidence_Interval`[i] <- paste(sort(D2_efficacies, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                Dose_Comparisons$`Mean_Combo_Efficacy_95%_Confidence_Interval`[i] <- paste(sort(Combo_efficacies, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
              #Calculating simulated IDAcomboscores and HRs
                if(Calculate_IDAcomboscore_And_Hazard_Ratio == TRUE){
                  #HRs
                    HRs_vs_D1 <- (1 - Combo_efficacies) / (1 - D1_efficacies)
                    HRs_vs_D2 <- (1 - Combo_efficacies) / (1 - D2_efficacies)
                  #HR standard errors
                    Dose_Comparisons$HR_vs_Control_Treatment_SE[i] <- sd(HRs_vs_D1)
                    Dose_Comparisons$HR_vs_Drug_to_Add_SE[i] <- sd(HRs_vs_D2)
                  #HR 95% confidence intervals
                    Dose_Comparisons$`HR_vs_Control_Treatment_95%_Confidence_Interval`[i] <- paste(sort(HRs_vs_D1, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                    Dose_Comparisons$`HR_vs_Drug_to_Add_95%_Confidence_Interval`[i] <- paste(sort(HRs_vs_D2, decreasing = FALSE)[c(index_2.5, index_97.5)], collapse = "_")
                  #Calculating HR p values
                    #Control therapy
                      #Calculating one-sided p-value with null hypothesis that HR >= 1
                        Dose_Comparisons$`p_HR_vs_Control_Treatment>=1`[i] <- sum(HRs_vs_D1 >= 1) / length(HRs_vs_D1)
                        #If p-value is 0, setting as p < minimum p value that can be estimated using this many simulations
                          if(Dose_Comparisons$`p_HR_vs_Control_Treatment>=1`[i] == 0){Dose_Comparisons$`p_HR_vs_Control_Treatment>=1`[i] <- paste0("<", 1/n_Simulations)}
                    #Drug_to_Add
                      #Calculating one-sided p-value with null hypothesis that HR >= 1
                        Dose_Comparisons$`p_HR_vs_Drug_to_Add>=1`[i] <- sum(HRs_vs_D2 >= 1) / length(HRs_vs_D2)
                        #If p-value is 0, setting as p < minimum p value that can be estimated using this many simulations
                          if(Dose_Comparisons$`p_HR_vs_Drug_to_Add>=1`[i] == 0){Dose_Comparisons$`p_HR_vs_Drug_to_Add>=1`[i] <- paste0("<", 1/n_Simulations)}
                  #IDAcomboscores
                    delta_viabilities <- pmin((1 - D1_efficacies), (1 - D2_efficacies)) - (1 - Combo_efficacies)
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
              rm(delta_viabilities, HR_C_over_Mbests, MC_IDAcomboscores, HRs_vs_D1, HRs_vs_D2)
            }
        }
    }

  #Returning Outputs
    #If Calculate_Uncertainty == FALSE, removing SE columns
      if(Calculate_Uncertainty == FALSE){
        Dose_Comparisons <- Dose_Comparisons[,-which(colnames(Dose_Comparisons) %in% c("Mean_Control_Treatment_Efficacy_SE", "Mean_Control_Treatment_Efficacy_95%_Confidence_Interval", "Mean_Drug_to_Add_Efficacy_SE", "Mean_Drug_to_Add_Efficacy_95%_Confidence_Interval", "Mean_Combo_Efficacy_SE", "Mean_Combo_Efficacy_95%_Confidence_Interval", "HR_vs_Control_Treatment_SE", "HR_vs_Control_Treatment_95%_Confidence_Interval", "p_HR_vs_Control_Treatment>=1", "HR_vs_Drug_to_Add_SE", "HR_vs_Drug_to_Add_95%_Confidence_Interval", "p_HR_vs_Drug_to_Add>=1", "IDA_Comboscore_SE", "IDA_Comboscore_95%_Confidence_Interval", "p_IDA_Comboscore<=0"))]
      }
    #If Calculate_IDAcomboscore_And_Hazard_Ratio == FALSE, removing HR and IDAcomboscore columns
      if(Calculate_IDAcomboscore_And_Hazard_Ratio == FALSE){
        Dose_Comparisons <- Dose_Comparisons[,-which(colnames(Dose_Comparisons) %in% c("HR_vs_Control_Treatment", "HR_vs_Control_Treatment_SE", "HR_vs_Control_Treatment_95%_Confidence_Interval", "p_HR_vs_Control_Treatment>=1", "HR_vs_Drug_to_Add", "HR_vs_Drug_to_Add_SE", "HR_vs_Drug_to_Add_95%_Confidence_Interval", "p_HR_vs_Drug_to_Add>=1", "IDA_Comboscore", "IDA_Comboscore_SE", "IDA_Comboscore_95%_Confidence_Interval", "p_IDA_Comboscore<=0"))]
      }
    #Replacing "Efficacy" with Efficacy_Metric_Name in column names of Dose_Comparisons
      colnames(Dose_Comparisons) <- gsub("Efficacy", Efficacy_Metric_Name, colnames(Dose_Comparisons))
    #Constructing Return_Object
      if(Return_Bootstrap_Values == FALSE){
        Return_Object <- list(Dose_Comparisons, gsub("^combo:", "", Control_Treatment_Name), Drug_to_Add, Usable_CellLines)
        names(Return_Object) <- c("Efficacy_Predictions", "Control_Treatment", "Drug_to_Add", "Cell_Lines_Used")
      } else if(Return_Bootstrap_Values == TRUE & Calculate_Uncertainty == TRUE){
        Return_Object <- list(Dose_Comparisons, gsub("^combo:", "", Control_Treatment_Name), Drug_to_Add, Usable_CellLines, D1_MC_Efficacies, D2_MC_Efficacies)
        names(Return_Object) <- c("Efficacy_Predictions", "Control_Treatment", "Drug_to_Add", "Cell_Lines_Used", "Bootstrap_Control_Efficacies", "Boostrap_Drug_to_Add_Efficacies")
      }
    #Returning output
      return(Return_Object)
}
