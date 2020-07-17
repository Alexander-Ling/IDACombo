#' Predicts and compares Bliss Independence efficacies for a pair of control and test treatments
#'
#' This function creates efficacy predictions for a pair of control and test treatments, each treatment consisting of a combination of one or more drugs, using monotherapy efficacy data and the assumptions of Bliss Independence. Concentrations must be specified for each drug in each treatment. IMPORTANT NOTE: This function is only applicable to drug efficacies measured on a scale from 0 to 1.
#'
#' @importFrom stats complete.cases rnorm sd aggregate
#'
#' @param Monotherapy_Data A data frame where each row contains information about the response of a single cell line to a single drug at a single concentration. Must minimally include columns containing the following information: cell line name, drug name, drug concentration, and measured drug efficacy. May optionally include a column recording the standard error (SE) of the measured drug efficacy.
#' @param Cell_Line_Name_Column A character vector of length 1 containing the name of the column in the Monotherapy_Data data frame which contains cell line names.
#' @param Drug_Name_Column A character vector of length 1 containing the name of the column in the Monotherapy_Data data frame which contains drug names.
#' @param Drug_Concentration_Column A character vector of length 1 containing the name of the column in the Monotherapy_Data data frame which contains drug concentrations.
#' @param Efficacy_Column A character vector of length 1 containing the name of the column in the Monotherapy_Data data frame which contains measured drug efficacies. Note, for Bliss Independence, efficacy must be expressed as a probability between 0 and 1.
#' @param LowerEfficacyIsBetterDrugEffect A logic vector of length 1 indicating whether or not lower values in Efficacy_Column indicate a more effective drug effect (i.e. for viability). Set TRUE if so. Otherwise, set FALSE if higher values in Efficacy_Column indicate a more effective drug response (i.e. for reduced viability).
#' @param EfficacyMetricName A character vector of length 1 indicating the name of the efficacy metric being used (i.e. Percent_Viability, Percent_Growth, etc.). Used to correctly label column names in output. Defaults to "Efficacy".
#' @param Control_Treatment_Drugs A character vector of length > 0 containing the names of the drugs in the control drug treatment for which efficacy predictions are to be made.
#' @param Control_Treatment_Drug_Concentrations A vector of drug concentrations for Control_Treatment_Drugs with the first concentration in Control_Treatment_Drug_Concentrations corresponding to the first drug in Control_Treatment_Drugs etc. Only one concentration may be specified for each drug in the control treatment, but, if a drug is included in both the control and test treatments, there is no need for the same concentration of that drug to be used in both treatments.
#' @param Test_Treatment_Drugs A character vector of length > 0 containing the names of the drugs in the control drug treatment for which efficacy predictions are to be made.
#' @param Test_Treatment_Drug_Concentrations A vector of drug concentrations for Test_Treatment_Drugs with the first concentration in Test_Treatment_Drug_Concentrations corresponding to the first drug in Test_Treatment_Drugs etc. Only one concentration may be specified for each drug in the test treatment, but, if a drug is included in both the control and test treatments, there is no need for the same concentration of that drug to be used in both treatments.
#' @param Calculate_Uncertainty A logic vector of length one indicating whether or not a Monte Carlo simulation should be performed to estimate uncertainties in the efficacy predictions based on uncertainties in the monotherapy efficacy measurements. Set TRUE if you wish to calculate uncertainties. Defaults to FALSE.
#' @param Efficacy_SE_Column A character vector of length 1 containing the name of the column in the Monotherapy_Data data frame which contains the standard errors of measured drug efficacies. Must be specified if Calculate_Uncertainty is set to TRUE.
#' @param n_Simulations A positive, non-zero integer vector of length 1 indicating the number of Monte Carlo iterations to be run when calculating output efficacy prediction uncertainties. Must be specified if Calculate_Uncertainty is set to TRUE. Defaults to 1000.
#' @param CalculateHazardRatio A logic vector of length 1 indicating whether or not a Hazard Ratio (HR) should be calculated between the control and test treatments. Set TRUE if so.
#' @param Average_Duplicate_Records A logic vector of length 1 indicating whether or not duplicated records (where a cell line has multiple records for being tested with a given drug at a given concentration) should be averaged. If TRUE, Efficacy values are averaged, and, if Calculate_Uncertainty is also TRUE, Efficacy_SE values are added in quadrature and divided by the number of duplicate records for that cell line/drug/concentration set.
#' @param Return_Monte_Carlo_HRs A logic vector of length 1 indicating whether or not the function should return the Hazard Ratios(HRs) simulated in a Monte Carlo simulation to estimate the Standard Erorr (SE) of the HR of the test treatment vs the control treatment. This parameter is ignored unless Calculate_Uncertainty = TRUE and CalculateHazardRatios = TRUE. This should be set to TRUE if the simulated HRs are needed to estimate uncertainties in downstream power analyses.
#'
#' @return If Return_Monte_Carlo_HRs = FALSE, this function returns a list with 4 elements: 1) A data frame with the produced efficacy predictions. 2) A data frame listing the control treatment drug names and concentrations. 3) A data frame listing the test treatment drug names and concentrations. 4) A character vector containing the names of the cell lines used to make the efficacy predictions. If Return_Monte_Carlo_HRs = TRUE & Calculate_Uncertainty = TRUE & CalculateHazardRatio = TRUE, this function returns a list with 5 elements: the first 4 elements being the same as when Return_Monte_Carlo_HRs = FALSE and the fifth element being a numeric vector of the Hazard Ratios simulated during the Monte Carlo simulation to estimate standard errors.
#'
#' @examples
#' #Loading Package
#'   library(IDACombo)
#'
#' #Making fake monotherapy dataset
#'   CellLineNames <- rep(c("CL1", "CL2", "CL3", "CL4", "CL5", "CL6"), 6)
#'   DrugNames <- c(rep("D1", 12), rep("D2", 12), rep("D3", 12))
#'   Concentrations <- c(rep(1, 6), rep(2, 6), rep(1.5, 6), rep(3, 6), rep("A", 6), rep("B", 6))
#'   Viability <- c(sample(seq(0.4,1,length.out = 10), 6, replace = TRUE),
#'                  sample(seq(0.2,0.8,length.out = 10), 6, replace = TRUE),
#'                  sample(seq(0.4,1,length.out = 10), 6, replace = TRUE),
#'                  sample(seq(0.2,0.6,length.out = 10), 6, replace = TRUE),
#'                  sample(seq(0.9,1,length.out = 10), 6, replace = TRUE),
#'                  sample(seq(0.2,0.6,length.out = 10), 6, replace = TRUE))
#'   Viability_SE <- Viability * sample(seq(0,0.1,length.out = 100), 36, replace = TRUE)
#'   Fake_Data <- data.frame(CellLineNames, DrugNames, Concentrations, Viability, Viability_SE)
#'
#' #Creating efficacy predictions for control and test treatments and comparing without
#' #uncertainty calculations
#'   #For case where drugs in test treatment are at reduced concentrations from
#'   #those used in the control treatment due to the addition of a third drug.
#'   #Note that this may mean that the test treatment is less effective than
#'   #the control treatment, such that the Hazard Ratio is > 1.
#'     BlissPredict.TestvsControl(Monotherapy_Data = Fake_Data,
#'                              Cell_Line_Name_Column = "CellLineNames",
#'                              Drug_Name_Column = "DrugNames",
#'                              Drug_Concentration_Column = "Concentrations",
#'                              Efficacy_Column = "Viability",
#'                              Control_Treatment_Drugs = c("D1", "D2"),
#'                              Control_Treatment_Drug_Concentrations = c(2, 3),
#'                              Test_Treatment_Drugs = c("D1", "D2", "D3"),
#'                              Test_Treatment_Drug_Concentrations = c(1, 1.5, "B"),
#'                              Calculate_Uncertainty = FALSE,
#'                              LowerEfficacyIsBetterDrugEffect = TRUE,
#'                              EfficacyMetricName = "Viability",
#'                              CalculateHazardRatio = TRUE,
#'                              Average_Duplicate_Records = FALSE)
#'
#'   #For case where drugs in test treatment are at same concentrations as
#'   #those used in the control treatment.
#'     BlissPredict.TestvsControl(Monotherapy_Data = Fake_Data,
#'                              Cell_Line_Name_Column = "CellLineNames",
#'                              Drug_Name_Column = "DrugNames",
#'                              Drug_Concentration_Column = "Concentrations",
#'                              Efficacy_Column = "Viability",
#'                              Control_Treatment_Drugs = c("D1", "D2"),
#'                              Control_Treatment_Drug_Concentrations = c(2, 3),
#'                              Test_Treatment_Drugs = c("D1", "D2", "D3"),
#'                              Test_Treatment_Drug_Concentrations = c(2, 3, "B"),
#'                              Calculate_Uncertainty = FALSE,
#'                              LowerEfficacyIsBetterDrugEffect = TRUE,
#'                              EfficacyMetricName = "Viability",
#'                              CalculateHazardRatio = TRUE,
#'                              Average_Duplicate_Records = FALSE)
#'
#' #Creating efficacy predictions for control and test treatments and comparing with
#' #uncertainty calculations but without returning Hazard Ratios that are generated
#' #using a Monte Carlo simulation to estimate standard errors.
#'
#'   BlissPredict.TestvsControl(Monotherapy_Data = Fake_Data,
#'                            Cell_Line_Name_Column = "CellLineNames",
#'                            Drug_Name_Column = "DrugNames",
#'                            Drug_Concentration_Column = "Concentrations",
#'                            Efficacy_Column = "Viability",
#'                            Control_Treatment_Drugs = c("D1", "D2"),
#'                            Control_Treatment_Drug_Concentrations = c(2, 3),
#'                            Test_Treatment_Drugs = c("D1", "D2", "D3"),
#'                            Test_Treatment_Drug_Concentrations = c(2, 3, "B"),
#'                            Calculate_Uncertainty = TRUE,
#'                            Efficacy_SE_Column = "Viability_SE",
#'                            n_Simulations = 1000,
#'                            LowerEfficacyIsBetterDrugEffect = TRUE,
#'                            EfficacyMetricName = "Viability",
#'                            CalculateHazardRatio = TRUE,
#'                            Average_Duplicate_Records = FALSE)
#'
#' #Creating efficacy predictions for control and test treatments and comparing with
#' #uncertainty calculations and with returning Hazard Ratios that are generated
#' #using a Monte Carlo simulation to estimate standard errors.
#'
#'   BlissPredict.TestvsControl(Monotherapy_Data = Fake_Data,
#'                            Cell_Line_Name_Column = "CellLineNames",
#'                            Drug_Name_Column = "DrugNames",
#'                            Drug_Concentration_Column = "Concentrations",
#'                            Efficacy_Column = "Viability",
#'                            Control_Treatment_Drugs = c("D1", "D2"),
#'                            Control_Treatment_Drug_Concentrations = c(2, 3),
#'                            Test_Treatment_Drugs = c("D1", "D2", "D3"),
#'                            Test_Treatment_Drug_Concentrations = c(2, 3, "B"),
#'                            Calculate_Uncertainty = TRUE,
#'                            Efficacy_SE_Column = "Viability_SE",
#'                            n_Simulations = 1000,
#'                            LowerEfficacyIsBetterDrugEffect = TRUE,
#'                            EfficacyMetricName = "Viability",
#'                            CalculateHazardRatio = TRUE,
#'                            Average_Duplicate_Records = FALSE,
#'                            Return_Monte_Carlo_HRs = TRUE)
#'
#' #Converting Viabilty to reduction in viability and redoing calculations
#' #without returning Monte Carlo Hazard Ratios. Note the change in the
#' #LowerEfficacyIsBetterDrugEffect flag from TRUE to FALSE
#'   Reduction_in_Viability <- 1-Viability
#'   Reduction_in_Viability_SE <- Viability_SE
#'   Fake_Data <- data.frame(CellLineNames,
#'                           DrugNames,
#'                           Concentrations,
#'                           Reduction_in_Viability,
#'                           Reduction_in_Viability_SE)
#'   BlissPredict.TestvsControl(Monotherapy_Data = Fake_Data,
#'                            Cell_Line_Name_Column = "CellLineNames",
#'                            Drug_Name_Column = "DrugNames",
#'                            Drug_Concentration_Column = "Concentrations",
#'                            Efficacy_Column = "Reduction_in_Viability",
#'                            Control_Treatment_Drugs = c("D1", "D2"),
#'                            Control_Treatment_Drug_Concentrations = c(2, 3),
#'                            Test_Treatment_Drugs = c("D1", "D2", "D3"),
#'                            Test_Treatment_Drug_Concentrations = c(2, 3, "B"),
#'                            Calculate_Uncertainty = TRUE,
#'                            Efficacy_SE_Column = "Reduction_in_Viability_SE",
#'                            n_Simulations = 1000,
#'                            LowerEfficacyIsBetterDrugEffect = FALSE,
#'                            EfficacyMetricName = "Reduction_in_Viability",
#'                            CalculateHazardRatio = TRUE,
#'                            Average_Duplicate_Records = FALSE)
#'
#' @export

BlissPredict.TestvsControl <- function(Monotherapy_Data, Cell_Line_Name_Column, Drug_Name_Column, Drug_Concentration_Column, Efficacy_Column, LowerEfficacyIsBetterDrugEffect, EfficacyMetricName = "Efficacy", Control_Treatment_Drugs, Control_Treatment_Drug_Concentrations, Test_Treatment_Drugs, Test_Treatment_Drug_Concentrations, Calculate_Uncertainty = FALSE, Efficacy_SE_Column = NULL, n_Simulations = 1000, CalculateHazardRatio = FALSE, Average_Duplicate_Records = FALSE, Return_Monte_Carlo_HRs = FALSE){
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
    if(! is.vector(EfficacyMetricName) | ! is.character(EfficacyMetricName) | ! length(EfficacyMetricName) == 1){
      stop("EfficacyMetricName is not a character vector of length 1.")
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
    if(! is.vector(Test_Treatment_Drugs) | ! is.character(Test_Treatment_Drugs) | ! length(Test_Treatment_Drugs) > 0){
      stop("Test_Treatment_Drugs is not a character vector with length > 0.")
    }
    if(! all(Test_Treatment_Drugs %in% Monotherapy_Data[,Drug_Name_Column])){
      missing.control.drugs <- Test_Treatment_Drugs[! Test_Treatment_Drugs %in% Monotherapy_Data[,Drug_Name_Column]]
      stop(paste0("No data for the following Test_Treatment_Drugs found in Monotherapy_Data: ", paste(missing.control.drugs, collapse = ", ")))
    }
    if(! is.vector(Test_Treatment_Drug_Concentrations) | ! length(Test_Treatment_Drug_Concentrations) == length(Test_Treatment_Drugs)){
      stop("Test_Treatment_Drug_Concentrations is not a vector with length = length(Test_Treatment_Drugs).")
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
    if(! is.vector(CalculateHazardRatio) | ! is.logical(CalculateHazardRatio) | ! length(CalculateHazardRatio) == 1){
      stop("CalculateHazardRatio is not a logical vector of length 1.")
    }
    if(! is.vector(Average_Duplicate_Records) | ! is.logical(Average_Duplicate_Records) | ! length(Average_Duplicate_Records) == 1){
      stop("Average_Duplicate_Records is not a logical vector of length 1.")
    }
    if(! is.vector(Return_Monte_Carlo_HRs) | ! is.logical(Return_Monte_Carlo_HRs) | ! length(Return_Monte_Carlo_HRs) == 1){
      stop("Return_Monte_Carlo_HRs is not a logical vector of length 1.")
    }

  #Organizing data into standard format based on column names provided for each desired set of information
  #Also subsetting to only include data pertaining to drugs in control or test treatment
    if(Calculate_Uncertainty == TRUE){
      Data <- Monotherapy_Data[Monotherapy_Data[,Drug_Name_Column] %in% c(Control_Treatment_Drugs, Test_Treatment_Drugs),c(Cell_Line_Name_Column, Drug_Name_Column, Drug_Concentration_Column, Efficacy_Column, Efficacy_SE_Column)]
      colnames(Data) <- c("CellLine", "Drug", "Conc", "Efficacy", "Efficacy_SE")
    } else {
      Data <- Monotherapy_Data[Monotherapy_Data[,Drug_Name_Column] %in% c(Control_Treatment_Drugs, Test_Treatment_Drugs),c(Cell_Line_Name_Column, Drug_Name_Column, Drug_Concentration_Column, Efficacy_Column)]
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

  #Throwing warning if Data$Efficacy has values < 0 or greater than 1
    if(min(Data$Efficacy) < 0 | max(Data$Efficacy) > 1){
      max.eff <- signif(max(Data$Efficacy), 3)
      min.eff <- signif(min(Data$Efficacy), 3)
      warning(paste0("Bliss Independence is only defined for efficacy as a probability from 0 to 1. Efficacy values in Monotherapy_Data fall between ", min.eff, " and ", max.eff, ". Efficacy values below zero will be rounded up to 0 and values above 1 will be rounded down to 1. Efficacy_SE values, if used, will not be changed."))
      Data$Efficacy <- pmin(Data$Efficacy, 1)
      Data$Efficacy <- pmax(Data$Efficacy, 0)
    }

  #Subsetting into control and test group data with specified concentrations
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
    TestData <- list(NULL)
    for(i in 1:length(Test_Treatment_Drugs)){
      TestData[[i]] <- Data[Data$Drug %in% Test_Treatment_Drugs[i],]
      #Checking that all provided concentrations are available for their respective drugs
        if(! Test_Treatment_Drug_Concentrations[i] %in% TestData[[i]]$Conc){
          stop(paste0(Test_Treatment_Drug_Concentrations[i], " concentration is unavailable for ", Test_Treatment_Drugs[i], " in test treatment."))
        } else {
          TestData[[i]] <- TestData[[i]][TestData[[i]]$Conc %in% Test_Treatment_Drug_Concentrations[i],]
        }
    }
    names(TestData) <- Test_Treatment_Drugs
    rm(Data)

  #Finding cell line overlap between all drugs
    ControlCellLines <- list(NULL)
    for(i in 1:length(ControlData)){
      ControlCellLines[[i]] <- sort(unique(ControlData[[i]]$CellLine))
    }
    TestCellLines <- list(NULL)
    for(i in 1:length(TestData)){
      TestCellLines[[i]] <- sort(unique(TestData[[i]]$CellLine))
    }
    All_Drug_CellLines <- c(ControlCellLines, TestCellLines)
    Usable_CellLines <- sort(unique(unlist(All_Drug_CellLines)))
    for(i in 1:length(All_Drug_CellLines)){
      Usable_CellLines <- Usable_CellLines[Usable_CellLines %in% All_Drug_CellLines[[i]]]
    }
    rm(ControlCellLines, TestCellLines, All_Drug_CellLines)

  #Checking again if at least 2 cell lines remain for all drugs. If not, exiting with NA
  #predictions and a warning.
    if(! length(Usable_CellLines) >= 2){
      #Returning NA predictions with warning due to too few cell lines.
      warning(paste0("<2 overlapping cell lines available for comparison of (", paste(Control_Treatment_Drugs, Control_Treatment_Drug_Concentrations, sep = "_", collapse = " + "), ") vs. (", paste(Test_Treatment_Drugs, Test_Treatment_Drug_Concentrations, sep = "_", collapse = " + "), ")"))
      if(Return_Monte_Carlo_HRs == FALSE){
        Return_Object <- list("Less than 2 overlapping cell lines available.", as.data.frame(cbind(Control_Treatment_Drugs, Control_Treatment_Drug_Concentrations), stringsAsFactors = F), as.data.frame(cbind(Test_Treatment_Drugs, Test_Treatment_Drug_Concentrations), stringsAsFactors = F), Usable_CellLines)
        names(Return_Object) <- c("Efficacy_Predictions", "Control_Treatment", "Test_Treatment", "Cell_Lines_Used")
        return(Return_Object)
      } else if(Return_Monte_Carlo_HRs == TRUE & CalculateHazardRatio == TRUE & Calculate_Uncertainty == TRUE){
        Return_Object <- list("Less than 2 overlapping cell lines available.", as.data.frame(cbind(Control_Treatment_Drugs, Control_Treatment_Drug_Concentrations), stringsAsFactors = F), as.data.frame(cbind(Test_Treatment_Drugs, Test_Treatment_Drug_Concentrations), stringsAsFactors = F), Usable_CellLines, "Less than 2 overlapping cell lines available.")
        names(Return_Object) <- c("Efficacy_Predictions", "Control_Treatment", "Test_Treatment", "Cell_Lines_Used", "Monte_Carlo_HRs")
        return(Return_Object)
      }
    }

  #Subsetting drug data to only include overlapping cell lines
    for(i in 1:length(ControlData)){
      ControlData[[i]] <- ControlData[[i]][ControlData[[i]]$CellLine %in% Usable_CellLines,]
    }
    for(i in 1:length(TestData)){
      TestData[[i]] <- TestData[[i]][TestData[[i]]$CellLine %in% Usable_CellLines,]
    }

  #Checking that cell lines aren't duplicated in each drug dataset
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
    for(i in 1:length(Test_Treatment_Drugs)){
      CL_Conc <- paste(TestData[[i]]$CellLine, TestData[[i]]$Conc, sep = "_")
      Dups <- CL_Conc[duplicated(CL_Conc)]
      if(length(Dups) > 0 & Average_Duplicate_Records == FALSE){
        warning(paste0("Duplicated information found for the following cell lines and ", Test_Treatment_Drugs[i], " concentrations in the test treatment. Average_Duplicate_Records = FALSE so duplicates removed: ", paste(Dups, collapse = ", ")))
        TestData[[i]] <- TestData[[i]][! duplicated(CL_Conc),]
      } else if(length(Dups) > 0 & Average_Duplicate_Records == TRUE){
        if(Calculate_Uncertainty == FALSE){
          #Simply averaging efficacy values
            colnames <- colnames(TestData[[i]])
            TestData[[i]] <- aggregate(TestData[[i]]$Efficacy, by = list(TestData[[i]]$CellLine, TestData[[i]]$Drug, TestData[[i]]$Conc), FUN = mean)
            colnames(TestData[[i]]) <- colnames
        } else if(Calculate_Uncertainty == TRUE){
          #Averaging efficacy
            Efficacy_Average <- aggregate(TestData[[i]]$Efficacy, by = list(TestData[[i]]$CellLine, TestData[[i]]$Drug, TestData[[i]]$Conc), FUN = mean)
            colnames(Efficacy_Average) <- c("CellLine", "Drug", "Conc", "Efficacy")
          #Calculating uncertainty in averaged efficacy by adding efficacy uncertainties in quadrature and dividing by number of values used in average
            Efficacy_Average_Uncertainties <- aggregate(TestData[[i]]$Efficacy_SE, by = list(TestData[[i]]$CellLine, TestData[[i]]$Drug, TestData[[i]]$Conc), FUN = function(x){sqrt(sum(x^2))/length(x)})
            colnames(Efficacy_Average_Uncertainties) <- c("CellLine", "Drug", "Conc", "Efficacy_SE")
          #Combining results
            TestData[[i]] <- merge(Efficacy_Average, Efficacy_Average_Uncertainties)
        }
      }
    }

  #Ordering cell lines the same for each drug
    for(i in 1:length(ControlData)){
      ControlData[[i]] <- ControlData[[i]][order(ControlData[[i]]$Conc, ControlData[[i]]$CellLine),]
    }
    for(i in 1:length(TestData)){
      TestData[[i]] <- TestData[[i]][order(TestData[[i]]$Conc, TestData[[i]]$CellLine),]
    }

  #Writing simple function to calculate expected efficacies based on Bliss Independence
  #Data should be a list of numeric efficacies
    Calc_Bliss_Independence_LowerEfficacyIsBetterDrugEffect <- function(Data){
      if(length(Data) > 1){
        for(i in 1:(length(Data)-1)){
          if(i == 1){
            Effect <- Data[[i]]*Data[[i+1]]
          } else {
            Effect <- Effect*Data[[i+1]]
          }
        }
      } else {
        Effect <- Data[[1]]
      }
      return(Effect)
    }
    Calc_Bliss_Independence_HigherEfficacyIsBetterDrugEffect <- function(Data){
      for(i in 1:length(Data)){
        Data[[i]] <- 1 - Data[[i]]
      }
      if(length(Data) > 1){
        for(i in 1:(length(Data)-1)){
          if(i == 1){
            Effect <- Data[[i]]*Data[[i+1]]
          } else {
            Effect <- Effect*Data[[i+1]]
          }
        }
      } else {
        Effect <- Data[[1]]
      }
      return(1 - Effect)
    }


  #Creating object to store prediction results in
    Prediction_Results <- as.data.frame(matrix(rep(NA, 6), nrow = 1))
    colnames(Prediction_Results) <- c("Mean_Control_Treatment_Efficacy", "Mean_Control_Treatment_Efficacy_SE", "Mean_Test_Treatment_Efficacy", "Mean_Test_Treatment_Efficacy_SE", "HR_Test_vs_Control_Treatment", "HR_Test_vs_Control_Treatment_SE")

  #Performing combination efficacy predictions for cases where lower efficacy values
  #indicate a more effective drug effect (i.e. efficacy = percent viability, percent growth, etc.)
    if(LowerEfficacyIsBetterDrugEffect == TRUE){
      #Calculating expected combination efficacy for each cell line
      #using Bliss Independence for both Control and Test treatment
        Control_Efficacy <- Calc_Bliss_Independence_LowerEfficacyIsBetterDrugEffect(lapply(ControlData, function(x){return(x$Efficacy)}))
        Test_Efficacy <- Calc_Bliss_Independence_LowerEfficacyIsBetterDrugEffect(lapply(TestData, function(x){return(x$Efficacy)}))
      #Calculating average efficacy accross all cell lines
        Prediction_Results$Mean_Control_Treatment_Efficacy <- mean(Control_Efficacy)
        Prediction_Results$Mean_Test_Treatment_Efficacy <- mean(Test_Efficacy)
        if(CalculateHazardRatio == TRUE){
          Prediction_Results$HR_Test_vs_Control_Treatment <- Prediction_Results$Mean_Test_Treatment_Efficacy / Prediction_Results$Mean_Control_Treatment_Efficacy
        }
        rm(Control_Efficacy, Test_Efficacy)

      #If Calculate_Uncertainty == TRUE, doing Monte Carlo simulation to estimate
      #uncertainties in output parameters based on uncertainties in monotherapy efficacies.
        if(Calculate_Uncertainty == TRUE){
          #Looping through each drug for each treatment and simulating efficacies based on
          #measured efficacies and SE's
            Control_MC_Efficacies <- as.list(NULL)
            for(i in 1:length(ControlData)){
              Control_MC_Efficacies[[i]] <- apply(ControlData[[i]][,c("Efficacy", "Efficacy_SE")], 1, function(x){rnorm(n = n_Simulations, mean = x[1], sd = x[2])})
              colnames(Control_MC_Efficacies[[i]]) <- ControlData[[i]]$CellLine
            }
            names(Control_MC_Efficacies) <- names(ControlData)
            Test_MC_Efficacies <- as.list(NULL)
            for(i in 1:length(TestData)){
              #Checking if this test drug is also in the control therapy.
              #If so, not randomly sampling for test viabilities. Instead, calculating how many SE's from
              #the measured value each simulated control viability fell, and matching that distance in
              #the simulated test viabilities. This is done because the drug in the control and test
              #therapies are not independent in this case--they will have come from the same dose-response curve.
                if(names(TestData)[i] %in% names(ControlData)){
                  #Calculating the number of SEs away from the measured value each simulated value is in the control therapy for this drug
                    Measured_Control_Data <- ControlData[[names(TestData)[i]]][,c("CellLine", "Efficacy", "Efficacy_SE")]
                    Measured_Control_Data <- Measured_Control_Data[match(colnames(Control_MC_Efficacies[[names(TestData)[i]]]), Measured_Control_Data$CellLine),]
                    Measured_Control_Efficacies <- matrix(Measured_Control_Data$Efficacy, ncol = length(Measured_Control_Data$Efficacy), nrow = n_Simulations, byrow = TRUE)
                    Measured_Control_Efficacy_SEs <- matrix(Measured_Control_Data$Efficacy_SE, ncol = length(Measured_Control_Data$Efficacy_SE), nrow = n_Simulations, byrow = TRUE)
                    SEs_deviated <- (Control_MC_Efficacies[[names(TestData)[i]]] - Measured_Control_Efficacies) / Measured_Control_Efficacy_SEs
                  #Calculating simulated test viabilities based on the SE deviations from the simulated control viabilities for this drug
                    Measured_Test_Data <- TestData[[names(TestData)[i]]][,c("CellLine", "Efficacy", "Efficacy_SE")]
                    Measured_Test_Data <- Measured_Test_Data[match(colnames(SEs_deviated), Measured_Test_Data$CellLine),]
                    Measured_Test_Efficacies <- matrix(Measured_Test_Data$Efficacy, ncol = length(Measured_Test_Data$Efficacy), nrow = n_Simulations, byrow = TRUE, dimnames = list(NULL, Measured_Test_Data$CellLine))
                    Measured_Test_Efficacy_SEs <- matrix(Measured_Test_Data$Efficacy_SE, ncol = length(Measured_Test_Data$Efficacy_SE), nrow = n_Simulations, byrow = TRUE, dimnames = list(NULL, Measured_Test_Data$CellLine))
                    Test_MC_Efficacies[[i]] <- Measured_Test_Efficacies + (SEs_deviated * Measured_Test_Efficacy_SEs)
                }
              #If test drug is not in control therapy, randomly sampling
                Test_MC_Efficacies[[i]] <- apply(TestData[[i]][,c("Efficacy", "Efficacy_SE")], 1, function(x){rnorm(n = n_Simulations, mean = x[1], sd = x[2])})
            }
            names(Test_MC_Efficacies) <- names(TestData)
            if(exists("Measured_Control_Data")){
              rm(Measured_Control_Data, Measured_Control_Efficacies, Measured_Control_Efficacy_SEs, Measured_Test_Data, Measured_Test_Efficacies, Measured_Test_Efficacy_SEs)
            }
          #Calculating expected combination efficacies for each cell line
          #using Bliss Independence for both Control and Test treatment
            Control_MC_Combo_Efficacies <- Calc_Bliss_Independence_LowerEfficacyIsBetterDrugEffect(Control_MC_Efficacies)
            Test_MC_Combo_Efficacies <- Calc_Bliss_Independence_LowerEfficacyIsBetterDrugEffect(Test_MC_Efficacies)
          #Calculating average efficacy accross all cell lines
            Mean_Control_Efficacies_MC <- rowMeans(Control_MC_Combo_Efficacies)
            Mean_Test_Efficacies_MC <- rowMeans(Test_MC_Combo_Efficacies)
            Prediction_Results$Mean_Control_Treatment_Efficacy_SE <- sd(Mean_Control_Efficacies_MC)
            Prediction_Results$Mean_Test_Treatment_Efficacy_SE <- sd(Mean_Test_Efficacies_MC)
            if(CalculateHazardRatio == TRUE){
              MC_HRs <- Mean_Test_Efficacies_MC / Mean_Control_Efficacies_MC
              Prediction_Results$HR_Test_vs_Control_Treatment_SE <- sd(MC_HRs)
            }
            rm(Control_MC_Combo_Efficacies, Test_MC_Combo_Efficacies, Mean_Test_Efficacies_MC, Mean_Control_Efficacies_MC)
        }
    }

  #Performing combination efficacy predictions for cases where lower efficacy values
  #indicate a less effective drug effect (i.e. efficacy = percent cell death, etc.)
    if(LowerEfficacyIsBetterDrugEffect == FALSE){
      #Calculating expected combination efficacy for each cell line
      #using Bliss Independence for both Control and Test treatment
        Control_Efficacy <- Calc_Bliss_Independence_HigherEfficacyIsBetterDrugEffect(lapply(ControlData, function(x){return(x$Efficacy)}))
        Test_Efficacy <- Calc_Bliss_Independence_HigherEfficacyIsBetterDrugEffect(lapply(TestData, function(x){return(x$Efficacy)}))
      #Calculating average efficacy accross all cell lines
        Prediction_Results$Mean_Control_Treatment_Efficacy <- mean(Control_Efficacy)
        Prediction_Results$Mean_Test_Treatment_Efficacy <- mean(Test_Efficacy)
        if(CalculateHazardRatio == TRUE){
          Prediction_Results$HR_Test_vs_Control_Treatment <- (1 - Prediction_Results$Mean_Test_Treatment_Efficacy) / (1 - Prediction_Results$Mean_Control_Treatment_Efficacy)
        }
        rm(Control_Efficacy, Test_Efficacy)

      #If Calculate_Uncertainty == TRUE, doing Monte Carlo simulation to estimate
      #uncertainties in output parameters based on uncertainties in monotherapy efficacies.
        if(Calculate_Uncertainty == TRUE){
          #Looping through each drug for each treatment and simulating efficacies based on
          #measured efficacies and SE's
            Control_MC_Efficacies <- as.list(NULL)
            for(i in 1:length(ControlData)){
              Control_MC_Efficacies[[i]] <- apply(ControlData[[i]][,c("Efficacy", "Efficacy_SE")], 1, function(x){rnorm(n = n_Simulations, mean = x[1], sd = x[2])})
              colnames(Control_MC_Efficacies[[i]]) <- ControlData[[i]]$CellLine
            }
            names(Control_MC_Efficacies) <- names(ControlData)
            Test_MC_Efficacies <- as.list(NULL)
            for(i in 1:length(TestData)){
              #Checking if this test drug is also in the control therapy.
              #If so, not randomly sampling for test viabilities. Instead, calculating how many SE's from
              #the measured value each simulated control viability fell, and matching that distance in
              #the simulated test viabilities. This is done because the drug in the control and test
              #therapies are not independent in this case--they will have come from the same dose-response curve.
                if(names(TestData)[i] %in% names(ControlData)){
                  #Calculating the number of SEs away from the measured value each simulated value is in the control therapy for this drug
                    Measured_Control_Data <- ControlData[[names(TestData)[i]]][,c("CellLine", "Efficacy", "Efficacy_SE")]
                    Measured_Control_Data <- Measured_Control_Data[match(colnames(Control_MC_Efficacies[[names(TestData)[i]]]), Measured_Control_Data$CellLine),]
                    Measured_Control_Efficacies <- matrix(Measured_Control_Data$Efficacy, ncol = length(Measured_Control_Data$Efficacy), nrow = n_Simulations, byrow = TRUE)
                    Measured_Control_Efficacy_SEs <- matrix(Measured_Control_Data$Efficacy_SE, ncol = length(Measured_Control_Data$Efficacy_SE), nrow = n_Simulations, byrow = TRUE)
                    SEs_deviated <- (Control_MC_Efficacies[[names(TestData)[i]]] - Measured_Control_Efficacies) / Measured_Control_Efficacy_SEs
                  #Calculating simulated test viabilities based on the SE deviations from the simulated control viabilities for this drug
                    Measured_Test_Data <- TestData[[names(TestData)[i]]][,c("CellLine", "Efficacy", "Efficacy_SE")]
                    Measured_Test_Data <- Measured_Test_Data[match(colnames(SEs_deviated), Measured_Test_Data$CellLine),]
                    Measured_Test_Efficacies <- matrix(Measured_Test_Data$Efficacy, ncol = length(Measured_Test_Data$Efficacy), nrow = n_Simulations, byrow = TRUE, dimnames = list(NULL, Measured_Test_Data$CellLine))
                    Measured_Test_Efficacy_SEs <- matrix(Measured_Test_Data$Efficacy_SE, ncol = length(Measured_Test_Data$Efficacy_SE), nrow = n_Simulations, byrow = TRUE, dimnames = list(NULL, Measured_Test_Data$CellLine))
                    Test_MC_Efficacies[[i]] <- Measured_Test_Efficacies + (SEs_deviated * Measured_Test_Efficacy_SEs)
                }
              #If test drug is not in control therapy, randomly sampling
                Test_MC_Efficacies[[i]] <- apply(TestData[[i]][,c("Efficacy", "Efficacy_SE")], 1, function(x){rnorm(n = n_Simulations, mean = x[1], sd = x[2])})
            }
            names(Test_MC_Efficacies) <- names(TestData)
            if(exists("Measured_Control_Data")){
              rm(Measured_Control_Data, Measured_Control_Efficacies, Measured_Control_Efficacy_SEs, Measured_Test_Data, Measured_Test_Efficacies, Measured_Test_Efficacy_SEs)
            }
          #Calculating expected combination efficacies for each cell line
          #using Bliss Independence for both Control and Test treatment
            Control_MC_Combo_Efficacies <- Calc_Bliss_Independence_HigherEfficacyIsBetterDrugEffect(Control_MC_Efficacies)
            Test_MC_Combo_Efficacies <- Calc_Bliss_Independence_HigherEfficacyIsBetterDrugEffect(Test_MC_Efficacies)
          #Calculating average efficacy accross all cell lines
            Mean_Control_Efficacies_MC <- rowMeans(Control_MC_Combo_Efficacies)
            Mean_Test_Efficacies_MC <- rowMeans(Test_MC_Combo_Efficacies)
            Prediction_Results$Mean_Control_Treatment_Efficacy_SE <- sd(Mean_Control_Efficacies_MC)
            Prediction_Results$Mean_Test_Treatment_Efficacy_SE <- sd(Mean_Test_Efficacies_MC)
            if(CalculateHazardRatio == TRUE){
              MC_HRs <- (1 - Mean_Test_Efficacies_MC) / (1 - Mean_Control_Efficacies_MC)
              Prediction_Results$HR_Test_vs_Control_Treatment_SE <- sd(MC_HRs)
            }
            rm(Control_MC_Combo_Efficacies, Test_MC_Combo_Efficacies, Mean_Test_Efficacies_MC, Mean_Control_Efficacies_MC)
        }
    }

  #Returning Outputs
    #If Calculate_Uncertainty == FALSE, removing SE columns
      if(Calculate_Uncertainty == FALSE){
        Prediction_Results <- Prediction_Results[,-which(colnames(Prediction_Results) %in% c("Mean_Control_Treatment_Efficacy_SE", "Mean_Test_Treatment_Efficacy_SE", "HR_Test_vs_Control_Treatment_SE"))]
      }
    #If CalculateHazardRatio == FALSE, removing HR columns
      if(CalculateHazardRatio == FALSE){
        Prediction_Results <- Prediction_Results[,-which(colnames(Prediction_Results) %in% c("HR_Test_vs_Control_Treatment", "HR_Test_vs_Control_Treatment_SE"))]
      }
    #Replacing "Efficacy" with EfficacyMetricName in column names of Dose_Comparisons
      colnames(Prediction_Results) <- gsub("Efficacy", EfficacyMetricName, colnames(Prediction_Results))
    #Constructing Return_Object
      if(Return_Monte_Carlo_HRs == FALSE){
        Return_Object <- list(Prediction_Results, as.data.frame(cbind(Control_Treatment_Drugs, Control_Treatment_Drug_Concentrations), stringsAsFactors = F), as.data.frame(cbind(Test_Treatment_Drugs, Test_Treatment_Drug_Concentrations), stringsAsFactors = F), Usable_CellLines)
        names(Return_Object) <- c("Efficacy_Predictions", "Control_Treatment", "Test_Treatment", "Cell_Lines_Used")
        return(Return_Object)
      } else if(Return_Monte_Carlo_HRs == TRUE & CalculateHazardRatio == TRUE & Calculate_Uncertainty == TRUE){
        Return_Object <- list(Prediction_Results, as.data.frame(cbind(Control_Treatment_Drugs, Control_Treatment_Drug_Concentrations), stringsAsFactors = F), as.data.frame(cbind(Test_Treatment_Drugs, Test_Treatment_Drug_Concentrations), stringsAsFactors = F), Usable_CellLines, MC_HRs)
        names(Return_Object) <- c("Efficacy_Predictions", "Control_Treatment", "Test_Treatment", "Cell_Lines_Used", "Monte_Carlo_HRs")
        return(Return_Object)
      }
    #Returning output
      return(Return_Object)
}
