#' @title pwrEWAS - A computationally efficient tool for compre-hensive power estimation in EWAS
#'
#' @description pwrEWAS is a computationally efficient tool to estimate power in EWAS as a function of sample and effect size for two-group comparisons of DNAm (e.g., case vs control, exposed vs non-exposed, etc.). Detailed description of in-/outputs, instructions and an example, as well as interpretations of the example results are provided in the following vignette: https://github.com/stefangraw/pwrEWAS/blob/master/vignettes/vignette.pdf
#' 
#' @param minTotSampleSize Minimum total sample size.
#' @param maxTotSampleSize Maximum total sample size.
#' @param SampleSizeSteps Sample size increments.
#' @param NcntPer Percentage sample group 1 (control group) (NcntPer = 0.5 indicates a balanced design).
#' @param targetDelta Target maximum difference in mean DNAm. (Either 'targetDelta' or 'deltaSD' should be specified)
#' @param deltaSD Standard deviation of simulated differences. (Either 'targetDelta' or 'deltaSD' should be specified)
#' @param J Number of CpGs tested/simulated (default: 100000).
#' @param targetDmCpGs Target number of DM CpGs.
#' @param tissueType Select a tissue type from the list of most commonly used tissue types: "Adult (PBMC)" (default), "Saliva", "Sperm", "Lymphoma", "Placenta", "Liver", "Colon", "Blood adult", "Blood 5 year olds", "Blood newborns", "Cord-blood (whole blood)" or "Cord-blood (PBMC)".
#' @param detectionLimit Smallest detectable difference in DNAm (default: 0.01).
#' @param DMmethod Method of Differential Methylation analysis: "limma" (default), "t-test (unequal var)", "t-test (equal var)", "Wilcox rank sum", "CPGassoc".
#' @param FDRcritVal FDRcritVal (default: 0.05).
#' @param core Number of threads for multi-threading (default: 1).
#' @param sims Number of simulated data sets (default: 50).
#'
#' @return pwrEWAS will return an object with the following four attributes: meanPower, powerArray, deltaArray, and metric, where metric contains marTypeI, classicalPower, FDR, and FDC
#'
#' @keywords DNAm microarray power
#' 
#' @export
#' 
#' @examples
#' outDelta = pwrEWAS(minTotSampleSize = 10,
#'                    maxTotSampleSize = 20,
#'                    SampleSizeSteps = 10,
#'                    NcntPer = 0.5,
#'                    targetDelta = c(0.2, 0.5),
#'                    J = 1000,
#'                    targetDmCpGs = 10,
#'                    tissueType = "Adult (PBMC)",
#'                    detectionLimit = 0.01,
#'                    DMmethod = "limma",
#'                    FDRcritVal = 0.05,
#'                    core = 2,
#'                    sims = 30)
#'                    
#' outSD = pwrEWAS(minTotSampleSize = 10,
#'                    maxTotSampleSize = 20,
#'                    SampleSizeSteps = 10,
#'                    NcntPer = 0.5,
#'                    deltaSD = c(0.02, 0.03),
#'                    J = 1000,
#'                    targetDmCpGs = 10,
#'                    tissueType = "Adult (PBMC)",
#'                    detectionLimit = 0.01,
#'                    DMmethod = "limma",
#'                    FDRcritVal = 0.05,
#'                    core = 2,
#'                    sims = 30)
pwrEWAS <- function(minTotSampleSize, # min total sample size
                    maxTotSampleSize, # max total sample size
                    SampleSizeSteps, # steps for increasing total sample size
                    NcntPer, # percentage of control sample
                    targetDelta = NULL, # vector of 99 percentile of the target max DM
                    deltaSD = NULL, # vector of sd(delta)
                    J = 100000, # number of simulated CpGs
                    targetDmCpGs, # target number for truely differentially methylated CpG
                    tissueType = c("Adult (PBMC)",
                                   "Saliva", 
                                   "Sperm", 
                                   "Lymphoma",
                                   "Placenta",
                                   "Liver",
                                   "Colon",
                                   "Blood adult",
                                   "Blood 5 year olds",
                                   "Blood newborns",
                                   "Cord-blood (whole blood)",
                                   "Cord-blood (PBMC)"), # tissue type that is used as reference to sample from
                    detectionLimit = 0.01, # lower bound for differential methylated CpG's to be considered TRULY differential methylated
                    DMmethod = c("limma", "t-test (unequal var)", "t-test (equal var)", "Wilcox rank sum", "CPGassoc"), # method to detect differential methylation
                    FDRcritVal = 0.05,
                    core = 1, # number of cores to multi thread
                    sims = 50
){
  # library(doParallel)
  # library(abind)
  # library(foreach)
  
  # # install non CRAN packages
  # if (!require("limma")){
  #   BiocManager::install("limma", version = "3.8", ask = FALSE)
  # }
  # if (!require("genefilter")){
  #   BiocManager::install("genefilter", version = "3.8", ask = FALSE)
  # }
  
  tissueType <- match.arg(tissueType)
  DMmethod <- match.arg(DMmethod)
  
  if(!is.null(targetDelta) & !is.null(deltaSD)) stop("Please specify only one: 'targetDelta' or 'deltaSD'")
  
  #load data
  methPara <- pwrEWAS.data:::loadDataset(tissueType)
  CpGonArray <- length(methPara$mu)
  
  output <- NULL
  totSampleSizes <- seq(minTotSampleSize, maxTotSampleSize, SampleSizeSteps)
  
  # initalize multi thread clusters
  cl <- parallel::makeCluster(core) # multi core
  doSNOW::registerDoSNOW(cl)
  # doParallel::registerDoParallel(cl) # multi core
  
  # combining function for foreach loops
  combine_tau <- function(listA, listB){
    if(is.null(listA)) return(listB)
    if(!is(listA[["power"]], "array") & !is(listA[["power"]], "matrix")) listA[["power"]] <- matrix(listA[["power"]])
    if(!is(listB[["power"]], "array") & !is(listB[["power"]], "matrix")) listB[["power"]] <- matrix(listB[["power"]])
    if(!is(listA[["metric"]]$marTypeI, "array") & !is(listA[["metric"]]$marTypeI, "matrix")) listA[["metric"]]$marTypeI <- matrix(listA[["metric"]]$marTypeI)
    if(!is(listB[["metric"]]$marTypeI, "array") & !is(listB[["metric"]]$marTypeI, "matrix")) listB[["metric"]]$marTypeI <- matrix(listB[["metric"]]$marTypeI)
    if(!is(listA[["metric"]]$classicalPower, "array") & !is(listA[["metric"]]$classicalPower, "matrix")) listA[["metric"]]$classicalPower <- matrix(listA[["metric"]]$classicalPower)
    if(!is(listB[["metric"]]$classicalPower, "array") & !is(listB[["metric"]]$classicalPower, "matrix")) listB[["metric"]]$classicalPower <- matrix(listB[["metric"]]$classicalPower)
    if(!is(listA[["metric"]]$FDR, "array") & !is(listA[["metric"]]$FDR, "matrix")) listA[["metric"]]$FDR <- matrix(listA[["metric"]]$FDR)
    if(!is(listB[["metric"]]$FDR, "array") & !is(listB[["metric"]]$FDR, "matrix")) listB[["metric"]]$FDR <- matrix(listB[["metric"]]$FDR)
    if(!is(listA[["metric"]]$FDC, "array") & !is(listA[["metric"]]$FDC, "matrix")) listA[["metric"]]$FDC <- matrix(listA[["metric"]]$FDC)
    if(!is(listB[["metric"]]$FDC, "array") & !is(listB[["metric"]]$FDC, "matrix")) listB[["metric"]]$FDC <- matrix(listB[["metric"]]$FDC)
    if(!is(listA[["metric"]]$probTP, "array") & !is(listA[["metric"]]$probTP, "matrix")) listA[["metric"]]$probTP <- matrix(listA[["metric"]]$probTP)
    if(!is(listB[["metric"]]$probTP, "array") & !is(listB[["metric"]]$probTP, "matrix")) listB[["metric"]]$probTP <- matrix(listB[["metric"]]$probTP)
    if(!is(listA[["delta"]], "list")) listA[["delta"]] <- list(listA[["delta"]])
    returnList <- list()
    returnList[["power"]] <- abind::abind(listA[["power"]], listB[["power"]], along = 3)
    returnList[["delta"]] <- listA[["delta"]]
    returnList[["delta"]][[length(listA[["delta"]])+1]] <- listB[["delta"]]
    returnList[["metric"]]$marTypeI <- abind::abind(listA[["metric"]]$marTypeI, listB[["metric"]]$marTypeI, along = 3)
    returnList[["metric"]]$classicalPower <- abind::abind(listA[["metric"]]$classicalPower, listB[["metric"]]$classicalPower, along = 3)
    returnList[["metric"]]$FDR <- abind::abind(listA[["metric"]]$FDR, listB[["metric"]]$FDR, along = 3)
    returnList[["metric"]]$FDC <- abind::abind(listA[["metric"]]$FDC, listB[["metric"]]$FDC, along = 3)
    returnList[["metric"]]$probTP <- abind::abind(listA[["metric"]]$probTP, listB[["metric"]]$probTP, along = 3)
    return(returnList)
  }
  
  combine_totSampleSizes <- function(listA, listB){
    if(is.null(listA)) return(listB)
    returnList <- list()
    returnList[["power"]] <- cbind(listA[["power"]], listB[["power"]])
    returnList[["delta"]] <- cbind(listA[["delta"]], listB[["delta"]])
    returnList[["metric"]]$marTypeI <- cbind(listA[["metric"]]$marTypeI, listB[["metric"]]$marTypeI)
    returnList[["metric"]]$classicalPower <- cbind(listA[["metric"]]$classicalPower, listB[["metric"]]$classicalPower)
    returnList[["metric"]]$FDR <- cbind(listA[["metric"]]$FDR, listB[["metric"]]$FDR)
    returnList[["metric"]]$FDC <- cbind(listA[["metric"]]$FDC, listB[["metric"]]$FDC)
    returnList[["metric"]]$probTP <- cbind(listA[["metric"]]$probTP, listB[["metric"]]$probTP)
    return(returnList)
  }
  
  
  if(is.null(deltaSD)){
    # finding tau and K for target DM CpG's
    cat(paste("[",Sys.time(),"] ", "Finding tau...", sep = ""))
    K <- NULL 
    tau <- NULL
    for(d in 1:length(targetDelta)){
      myTau <- getTau(targetDmCpGs, targetDelta[d], methPara, detectionLimit, J, CpGonArray)
      tau[d] <- myTau$tau
      K[d] <- myTau$K # number of changed CpGs
    }
    cat(paste("done", " [",Sys.time(),"]\n", sep = ""))
    print(paste("The following taus were chosen: ", paste(tau, collapse = ", "), sep = ""))
  } else {
    tau <- deltaSD 
    K <- NULL
    for(d in 1:length(tau)){
      K[d] <- getK(targetDmCpGs, methPara, detectionLimit, J, CpGonArray, tau)
    }
  }
  
  
  
  # main function
  startTime = Sys.time()
  cat(paste("[",startTime,"] ", "Running simulation\n", sep = ""))
  iterations <- length(totSampleSizes) * length(tau)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  multiThreadOut <- foreach(d = 1:length(tau), 
                            .combine = combine_tau,
                            .packages=c("truncnorm", "limma", "CpGassoc", "genefilter"),
                            .export = c("getAlphBet", "getMeanVar", "beta2Mvalue", "limma", "ttestSlow", "ttestFast", "Wilcox", "CPGassoc")) %:%
    foreach(Ntot = totSampleSizes, .combine = combine_totSampleSizes, .options.snow = opts) %dopar% { 
      
      setTxtProgressBar(pb, (d-1)*length(totSampleSizes) + which(Ntot==totSampleSizes))
      
      Ncnt <- round(Ntot * NcntPer)
      Ntx <- Ntot - Ncnt
      marPower <- NULL
      deltaSim <- NULL
      
      marTypeI <- NULL
      FDR <- NULL
      classicalPower <- NULL
      FDC <- NULL
      probTP <- NULL
      
      for(sim in 1:sims){
        
        ## sample CpGs
        cpgIdx <- sample(x = 1:CpGonArray, size = J, replace = TRUE) # pick J random CpG's to be simulated
        cpgIdxName <- paste(1:J, "_", rownames(methPara)[cpgIdx], sep = "") # ensuring unique CpG name (allowing unique sampling with replacement)
        changedCpgsIdx <- sample(x = cpgIdx, size = K[d]) # pick K random CpG's to be changed in mean meth
        changedCpgsIdxName <- cpgIdxName[match(changedCpgsIdx, cpgIdx)]
        
        ## Change Mu for "changedCpgsIdx"'s CpG's
        # drawing delta from truncated normal
        delta <- truncnorm::rtruncnorm(1, mean = 0, sd = as.numeric(tau[d]), 
                                       a=0.5 - methPara$mu[changedCpgsIdx] - sqrt(0.25-methPara$var[changedCpgsIdx]), 
                                       b=0.5 - methPara$mu[changedCpgsIdx] + sqrt(0.25-methPara$var[changedCpgsIdx]))
        deltaSim <- c(deltaSim, delta)  # multi core
        meaningfulDM <- (abs(delta) >= detectionLimit)
        meaningfulDMName <- changedCpgsIdxName[meaningfulDM]
        
        # changing mean methylation for specific CpG's (changedCpgsIdx)
        muToBeSimuUNchanged <- methPara$mu[cpgIdx]
        muToBeSimuChanged   <- methPara$mu[cpgIdx]
        muToBeSimuChanged[match(changedCpgsIdx, cpgIdx)] <- muToBeSimuChanged[match(changedCpgsIdx, cpgIdx)] + delta
        
        ## get alpha and beta 
        params_unchanged <- getAlphBet(myMean = muToBeSimuUNchanged, myVar = methPara$var[cpgIdx])
        alpha_unchanged <- params_unchanged$alpha
        beta_unchanged <- params_unchanged$beta
        params_changed <- getAlphBet(myMean = muToBeSimuChanged, myVar = methPara$var[cpgIdx])
        alpha_changed <- params_changed$alpha
        beta_changed <- params_changed$beta
        
        ## simulate baseline beta values
        g1Beta <- NULL
        g2Beta <- NULL
        g1Beta <- matrix(rbeta(J*Ncnt, rep(alpha_unchanged, each = Ncnt), rep(beta_unchanged, each = Ncnt)), ncol = Ncnt, byrow = TRUE) 
        g2Beta <- matrix(rbeta(J*Ntx, rep(alpha_changed, each = Ntx), rep(beta_changed, each = Ntx)), ncol = Ntx, byrow = TRUE) 
        g1Beta[g1Beta == 1] <- max(g1Beta[g1Beta != 1]) # replacing 0/1 by min/max
        g2Beta[g2Beta == 1] <- max(g2Beta[g2Beta != 1])
        g1Beta[g1Beta == 0] <- min(g1Beta[g1Beta != 0])
        g2Beta[g2Beta == 0] <- min(g2Beta[g2Beta != 0])
        rownames(g1Beta) <- rownames(g2Beta) <- paste(1:J,"_",names(alpha_unchanged),sep = "")
        
        # Tests 
        ## t-test slow (unequal var)
        if(DMmethod == "t-test (unequal var)") {
          DMtest <- ttestSlow(g1Beta,g2Beta,Ncnt,Ntx,paired=FALSE)
        } else 
          
          ## t-test (equal var)
          if(DMmethod == "t-test (equal var)") {
            DMtest <- ttestFast(g1Beta,g2Beta,Ncnt,Ntx)
          } else
            
            ## CPG assoc
            if(DMmethod == "CPGassoc") {
              DMtest <- CPGassoc(g1Beta, g2Beta, Ncnt,Ntx)
            } else
              
              ## Wilcox rank sum test
              if(DMmethod == "Wilcox rank sum") {
                DMtest <- Wilcox(g1Beta, g2Beta, Ncnt,Ntx)
              } else 
                
                ## limma
                if(DMmethod == "limma") {
                  DMtest <- limma(g1Beta, g2Beta, Ncnt,Ntx)
                }
        
        # table from paper
        notDM  <- cpgIdxName[!(cpgIdxName %in% changedCpgsIdxName)]
        DM_negligible <- changedCpgsIdxName[!(changedCpgsIdxName %in% meaningfulDMName)]
        DM_meaningful <- changedCpgsIdxName[changedCpgsIdxName %in% meaningfulDMName]
        
        FP  <- intersect(cpgIdxName[DMtest$fdr<FDRcritVal], notDM) # FP
        NP <- intersect(cpgIdxName[DMtest$fdr<FDRcritVal], DM_negligible) # NP
        TP <- intersect(cpgIdxName[DMtest$fdr<FDRcritVal], DM_meaningful) # TP
        detectedCpGs <- cpgIdxName[DMtest$fdr<FDRcritVal]
        
        TN <- intersect(cpgIdxName[!(DMtest$fdr<FDRcritVal)], notDM) # TN
        NN <- intersect(cpgIdxName[!(DMtest$fdr<FDRcritVal)], DM_negligible) # NN
        FN <- intersect(cpgIdxName[!(DMtest$fdr<FDRcritVal)], DM_meaningful) # FN
        
        
        ## metrics
        marPower[sim] <- ifelse(length(DM_meaningful) > 0, length(TP)/length(DM_meaningful), NA)
        marTypeI[sim] <- ifelse(length(notDM) > 0, length(FP)/length(notDM), NA) 
        FDR[sim] <- ifelse(length(detectedCpGs) > 0, length(FP)/length(detectedCpGs), NA)
        FDC[sim] <- ifelse(length(TP) > 0, (length(FP))/length(TP), NA) # expected false discovery for each true discovery
        classicalPower[sim] <- (length(NP)+length(TP))/(length(DM_negligible)+length(DM_meaningful))
        probTP[sim] <- ifelse(length(TP)>0, 1, 0)
        
      } # end sim
      
      outSim=list() 
      outSim[["power"]] <- marPower # multi core
      outSim[["delta"]] <- deltaSim # multi core
      outSim[["metric"]]$marTypeI <- marTypeI
      outSim[["metric"]]$FDR <- FDR
      outSim[["metric"]]$classicalPower <- classicalPower
      outSim[["metric"]]$FDC <- FDC
      outSim[["metric"]]$probTP <- probTP
      outSim
    } # end tau and totSampleSizes
  close(pb)
  parallel::stopCluster(cl) 
  cat(paste("[",startTime,"] Running simulation ... done [",Sys.time(),"]\n", sep = ""))
  
  # dirty fix
  if(is.null(targetDelta)) targetDelta <- tau
  
  # calculate mean marginal power and adding names 
  if(length(targetDelta) == 1 & length(totSampleSizes) == 1)  output$meanPower <- matrix(mean(multiThreadOut[["power"]], na.rm=TRUE))
  if(length(targetDelta) == 1 & length(totSampleSizes) > 1)   output$meanPower <- matrix(apply(multiThreadOut[["power"]], 2, mean, na.rm=TRUE), ncol = 1)
  if(length(targetDelta) > 1)                                 output$meanPower <- apply(multiThreadOut[["power"]], c(2,3), mean, na.rm=TRUE)
  rownames(output$meanPower) <- totSampleSizes
  colnames(output$meanPower) <- targetDelta
  
  # powerArray
  if(length(targetDelta) == 1) output$powerArray <- array(data = multiThreadOut[["power"]], dim = c(sims, length(totSampleSizes), length(targetDelta)))
  if(length(targetDelta) > 1 & length(totSampleSizes) == 1) output$powerArray <- multiThreadOut[["power"]]
  if(length(targetDelta) > 1) output$powerArray <- multiThreadOut[["power"]]
  dimnames(output$powerArray) <- list(1:sims, totSampleSizes, targetDelta)  
  
  # deltaArray
  if(length(targetDelta) == 1 & length(totSampleSizes) == 1)  output$deltaArray <- list(matrix(multiThreadOut[["delta"]]))
  if(length(targetDelta) == 1 & length(totSampleSizes) > 1)   output$deltaArray <- list(multiThreadOut[["delta"]])
  if(length(targetDelta) > 1 & length(totSampleSizes) == 1)   output$deltaArray <- lapply(multiThreadOut[["delta"]],as.matrix)
  if(length(targetDelta) > 1 & length(totSampleSizes) > 1)    output$deltaArray <- multiThreadOut[["delta"]]
  names(output$deltaArray) <- targetDelta
  for(d in 1:length(targetDelta)){
    colnames(output$deltaArray[[d]]) <- totSampleSizes
  }
  
  
  # calculate mean marTypeI and adding names 
  if(length(targetDelta) == 1 & length(totSampleSizes) == 1)  output$metric$marTypeI <- matrix(mean(multiThreadOut[["metric"]]$marTypeI, na.rm=TRUE))
  if(length(targetDelta) == 1 & length(totSampleSizes) > 1)   output$metric$marTypeI <- matrix(apply(multiThreadOut[["metric"]]$marTypeI, 2, mean, na.rm=TRUE), ncol = 1)
  if(length(targetDelta) > 1)                                 output$metric$marTypeI <- apply(multiThreadOut[["metric"]]$marTypeI, c(2,3), mean, na.rm=TRUE)
  rownames(output$metric$marTypeI) <- totSampleSizes
  colnames(output$metric$marTypeI) <- targetDelta
  
  # calculate mean classicalPower and adding names 
  if(length(targetDelta) == 1 & length(totSampleSizes) == 1)  output$metric$classicalPower <- matrix(mean(multiThreadOut[["metric"]]$classicalPower, na.rm=TRUE))
  if(length(targetDelta) == 1 & length(totSampleSizes) > 1)   output$metric$classicalPower <- matrix(apply(multiThreadOut[["metric"]]$classicalPower, 2, mean, na.rm=TRUE), ncol = 1)
  if(length(targetDelta) > 1)                                 output$metric$classicalPower <- apply(multiThreadOut[["metric"]]$classicalPower, c(2,3), mean, na.rm=TRUE)
  rownames(output$metric$classicalPower) <- totSampleSizes
  colnames(output$metric$classicalPower) <- targetDelta
  
  # calculate mean FDR and adding names 
  if(length(targetDelta) == 1 & length(totSampleSizes) == 1)  output$metric$FDR <- matrix(mean(multiThreadOut[["metric"]]$FDR, na.rm=TRUE))
  if(length(targetDelta) == 1 & length(totSampleSizes) > 1)   output$metric$FDR <- matrix(apply(multiThreadOut[["metric"]]$FDR, 2, mean, na.rm=TRUE), ncol = 1)
  if(length(targetDelta) > 1)                                 output$metric$FDR <- apply(multiThreadOut[["metric"]]$FDR, c(2,3), mean, na.rm=TRUE)
  rownames(output$metric$FDR) <- totSampleSizes
  colnames(output$metric$FDR) <- targetDelta
  
  # calculate mean FDC and adding names 
  if(length(targetDelta) == 1 & length(totSampleSizes) == 1)  output$metric$FDC <- matrix(mean(multiThreadOut[["metric"]]$FDC, na.rm=TRUE))
  if(length(targetDelta) == 1 & length(totSampleSizes) > 1)   output$metric$FDC <- matrix(apply(multiThreadOut[["metric"]]$FDC, 2, mean, na.rm=TRUE), ncol = 1)
  if(length(targetDelta) > 1)                                 output$metric$FDC <- apply(multiThreadOut[["metric"]]$FDC, c(2,3), mean, na.rm=TRUE)
  rownames(output$metric$FDC) <- totSampleSizes
  colnames(output$metric$FDC) <- targetDelta
  
  # calculate probability of detecting at least one TP and adding names 
  if(length(targetDelta) == 1 & length(totSampleSizes) == 1)  output$metric$probTP <- matrix(mean(multiThreadOut[["metric"]]$probTP, na.rm=TRUE))
  if(length(targetDelta) == 1 & length(totSampleSizes) > 1)   output$metric$probTP <- matrix(apply(multiThreadOut[["metric"]]$probTP, 2, mean, na.rm=TRUE), ncol = 1)
  if(length(targetDelta) > 1)                                 output$metric$probTP <- apply(multiThreadOut[["metric"]]$probTP, c(2,3), mean, na.rm=TRUE)
  rownames(output$metric$probTP) <- totSampleSizes
  colnames(output$metric$probTP) <- targetDelta
  
  return(output)
}
