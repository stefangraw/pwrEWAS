#' @title pwrEWAS - A computationally efficient tool for compre-hensive power estimation in EWAS
#'
#' @description pwrEWAS is a computationally efficient tool to estimate power in EWAS as a function of sample and effect size for two-group comparisons of DNAm (e.g., case vs control, exposed vs non-exposed, etc.). Detailed description of in-/outputs, instructions and an example, as well as interpretations of the example results are provided in the following vignette:
#' 
#' @param minTotSampleSize Minimum total sample size (default: 10).
#' @param maxTotSampleSize Maximum total sample size (default: 50).
#' @param SampleSizeSteps Sample size increments (default: 10).
#' @param NcntPer Percentage sample group 1 (control group) (default: 0.5).
#' @param targetDelta Target maximum difference in mean DNAm (default: c(0.2, 0.5).
#' @param J Number of CpGs tested/simulated (default: 10000).
#' @param targetDmCpGs Target number of DM CpGs (default: 100).
#' @param tissueType Select a tissue type from the list of most commonly used tissue types: "Adult (PBMC)" (default), "Saliva", "Lymphoma", "Placenta", "Liver", "Colon", "Peripheral Leukocytes", "Blood 5 year olds", "Blood newborns", "Cord-blood (whole blood)" or "Cord-blood (PBMC)".
#' @param detectionLimit Smallest detectable difference in DNAm (default: 0.01).
#' @param DMmethod Method of Differential Methylation analysis: "limma" (default), "t-test (unequal var)", "t-test (equal var)", "Wilcox rank sum", "CPGassoc".
#' @param FDRcritVal FDRcritVal (default: 0.05).
#' @param core Number of threads for multi-threading (default: 1).
#' @param sim Number of simulated data sets (default: 50).
#'
#' @return pwrEWAS will return an object with the following three attributes: meanPower, powerArray, deltaArray
#'
#' @keywords DNAm microarray power
#' 
#' @export
#' 
#' @examples
#' pwrEWAS(minTotSampleSize = 10,
#'                    maxTotSampleSize = 50,
#'                    SampleSizeSteps = 10,
#'                    NcntPer = 10,
#'                    targetDelta = c(0.2, 0.5),
#'                    J = 10000,
#'                    targetDmCpGs = 100,
#'                    tissueType = "Adult (PBMC)",
#'                    detectionLimit = 0.01,
#'                    DMmethod = "limma",
#'                    FDRcritVal = 0.05,
#'                    core = 4,
#'                    sims = 50)

pwrEWAS = function(minTotSampleSize = 10, # min total sample size
                   maxTotSampleSize = 50, # max total sample size
                   SampleSizeSteps = 10, # steps for increasing total sample size
                   NcntPer = 0.5, # percentage of control sample
                   targetDelta = c(0.2, 0.5), # vector of 99 percentile of the target max DM
                   # CpGonArray, # "max" CpG's on array
                   J = 10000, # number of simulated CpGs
                   targetDmCpGs = 100, # target number for truely differentially methylated CpG
                   # methPara, # dataframe with mean- and var-vector
                   tissueType = c("Adult (PBMC)",
                                  "Saliva", 
                                  "Lymphoma",
                                  "Placenta",
                                  "Liver",
                                  "Colon",
                                  "Peripheral Leukocytes",
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
  library(doParallel)
  library(abind)
  library(foreach)
  # source(file = "functions.R")  # multi core
  tissueType = match.arg(tissueType)
  DMmethod = match.arg(DMmethod)
  
  #load data
  methPara = loadDataset(tissueType)
  CpGonArray = length(methPara$mu)
  
  output = NULL
  totSampleSizes = seq(minTotSampleSize, maxTotSampleSize, SampleSizeSteps)
  
  # initalize multi thread clusters
  cl <- makeCluster(core) # multi core
  registerDoParallel(cl) # multi core
  
  
  # combining function for foreach loops
  combine_tau = function(listA, listB){
    if(class(listA[["power"]]) !="array" & class(listA[["power"]]) !="matrix") listA[["power"]] = matrix(listA[["power"]])
    if(class(listB[["power"]]) !="array" & class(listB[["power"]]) !="matrix") listB[["power"]] = matrix(listB[["power"]])
    if(class(listA[["delta"]]) !="list") listA[["delta"]] = list(listA[["delta"]])
    returnList = list()
    returnList[["power"]] = abind(listA[["power"]], listB[["power"]], along = 3)
    returnList[["delta"]] = listA[["delta"]]
    returnList[["delta"]][[length(listA[["delta"]])+1]] = listB[["delta"]]
    return(returnList)
  }
  combine_totSampleSizes = function(listA, listB){
    returnList = list()
    returnList[["power"]] = cbind(listA[["power"]], listB[["power"]])
    returnList[["delta"]] = cbind(listA[["delta"]], listB[["delta"]])
    return(returnList)
  }
  
  # finding tau and K for target DM CpG's
  cat(paste("[",Sys.time(),"] ", "Finding tau...", sep = ""))
  K = NULL; tau = NULL
  for(d in 1:length(targetDelta)){
    myTau = getTau(targetDmCpGs, targetDelta[d], methPara, detectionLimit, J, CpGonArray)
    tau[d] = myTau$tau
    K[d] = myTau$K # number of changed CpGs
  }
  cat(paste("done", " [",Sys.time(),"]\n", sep = ""))
  # cat("done\n")
  
  # main function
  cat(paste("[",Sys.time(),"] ", "Running simulation ...", sep = ""))
  multiThreadOut = foreach(d = 1:length(tau), 
                           .combine = combine_tau,
                           .packages=c("truncnorm", "limma", "CpGassoc", "genefilter")) %:%
    foreach(Ntot = totSampleSizes, .combine = combine_totSampleSizes) %dopar% { 
      
      # source(file = "functions.R")  # multi core
      Ncnt = round(Ntot * NcntPer)
      Ntx = Ntot - Ncnt
      powerSim = NULL
      deltaSim = NULL
      
      for(sim in 1:sims){
        
        ## sample CpGs
        cpgIdx = sample(x = 1:CpGonArray, size = J, replace = T) # pick J random CpG's to be simulated
        cpgIdxName = paste(1:J, "_", names(methPara$mu)[cpgIdx], sep = "") # ensuring unique CpG name (allowing unique sampling with replacement)
        changedCpgsIdx = sample(x = cpgIdx, size = K[d]) # pick K random CpG's to be changed in mean meth
        changedCpgsIdxName = cpgIdxName[match(changedCpgsIdx, cpgIdx)]
        
        ## Change Mu for "changedCpgsIdx"'s CpG's
        # drawing delta from truncated normal
        delta = rtruncnorm(1, mean = 0, sd = as.numeric(tau[d]), 
                           a=0.5 - methPara$mu[changedCpgsIdx] - sqrt(0.25-methPara$var[changedCpgsIdx]), 
                           b=0.5 - methPara$mu[changedCpgsIdx] + sqrt(0.25-methPara$var[changedCpgsIdx]))
        deltaSim = c(deltaSim, delta)  # multi core
        meaningfulDM = (abs(delta) >= detectionLimit)
        meaningfulDMName = changedCpgsIdxName[meaningfulDM]
        
        # changing mean methylation for specific CpG's (changedCpgsIdx)
        muToBeSimuUNchanged = methPara$mu[cpgIdx]
        muToBeSimuChanged   = methPara$mu[cpgIdx]
        muToBeSimuChanged[match(changedCpgsIdx, cpgIdx)] = muToBeSimuChanged[match(changedCpgsIdx, cpgIdx)] + delta
        
        ## get alpha and beta 
        params_unchanged = getAlphBet(myMean = muToBeSimuUNchanged, myVar = methPara$var[cpgIdx])
        alpha_unchanged = params_unchanged$alpha
        beta_unchanged = params_unchanged$beta
        params_changed = getAlphBet(myMean = muToBeSimuChanged, myVar = methPara$var[cpgIdx])
        alpha_changed = params_changed$alpha
        beta_changed = params_changed$beta
        
        ## simulate baseline beta values
        g1Beta = NULL
        g2Beta = NULL
        g1Beta = matrix(rbeta(J*Ncnt, rep(alpha_unchanged, each = Ncnt), rep(beta_unchanged, each = Ncnt)), ncol = Ncnt, byrow = TRUE) 
        g2Beta = matrix(rbeta(J*Ntx, rep(alpha_changed, each = Ntx), rep(beta_changed, each = Ntx)), ncol = Ntx, byrow = TRUE) 
        g1Beta[g1Beta == 1] = max(g1Beta[g1Beta != 1]) # replacing 0/1 by min/max
        g2Beta[g2Beta == 1] = max(g2Beta[g2Beta != 1])
        g1Beta[g1Beta == 0] = min(g1Beta[g1Beta != 0])
        g2Beta[g2Beta == 0] = min(g2Beta[g2Beta != 0])
        rownames(g1Beta) = rownames(g2Beta) = paste(1:J,"_",names(alpha_unchanged),sep = "")
        
        # Tests 
        ## t-test slow (unequal var)
        if(DMmethod == "t-test (unequal var)") {
          DMtest = ttestSlow(g1Beta,g2Beta,Ncnt,Ntx,paired=FALSE)
        } else 
          
          ## t-test (equal var)
          if(DMmethod == "t-test (equal var)") {
            DMtest = ttestFast(g1Beta,g2Beta,Ncnt,Ntx)
          } else
            
            ## CPG assoc
            if(DMmethod == "CPGassoc") {
              DMtest = CPGassoc(g1Beta, g2Beta, Ncnt,Ntx)
            } else
              
              ## Wilcox rank sum test
              if(DMmethod == "Wilcox rank sum") {
                DMtest = Wilcox(g1Beta, g2Beta, Ncnt,Ntx)
              } else 
                
                ## limma
                if(DMmethod == "limma") {
                  DMtest = limma(g1Beta, g2Beta, Ncnt,Ntx)
                }
        
        # table from paper
        G0  = cpgIdxName[!(cpgIdxName %in% changedCpgsIdxName)]
        G1A = changedCpgsIdxName[!(changedCpgsIdxName %in% meaningfulDMName)]
        G1B = changedCpgsIdxName[changedCpgsIdxName %in% meaningfulDMName]
        
        V  = intersect(cpgIdxName[DMtest$fdr<FDRcritVal], G0)
        Sa = intersect(cpgIdxName[DMtest$fdr<FDRcritVal], G1A)
        Sb = intersect(cpgIdxName[DMtest$fdr<FDRcritVal], G1B)
        R = cpgIdxName[DMtest$fdr<FDRcritVal]
        
        G0_V = intersect(cpgIdxName[!(DMtest$fdr<FDRcritVal)], G0)
        G1A_Sa = intersect(cpgIdxName[!(DMtest$fdr<FDRcritVal)], G1A)
        G1B_Sb = intersect(cpgIdxName[!(DMtest$fdr<FDRcritVal)], G1B)
        
        
        ## metrics
        if(length(G0) > 0) marTypeI = length(V)/length(G0) else marTypeI = 1
        if(length(R) > 0) FDR = length(V)/length(R) else FDR = 1
        power = (length(Sa)+length(Sb))/(length(G1A)+length(G1B))
        if(length(G1B) > 0) marPower = length(Sb)/length(G1B) else marPower = 0
        if(length(Sb) > 0) FDC = (length(V)+length(Sa))/length(Sb) else FDC = Inf # expected false discovery for each true discovery
        
        powerSim[sim] = marPower # multi core
      } # end sim
      
      outSim=list() 
      outSim[["power"]] = powerSim # multi core
      outSim[["delta"]] = deltaSim # multi core
      outSim
    } # end tau and totSampleSizes
  stopCluster(cl) 
  cat(paste("done", " [",Sys.time(),"]\n", sep = ""))
  
  # calculate mean and adding names 
  if(length(targetDelta) == 1 & length(totSampleSizes) == 1)  output$meanPower = matrix(mean(multiThreadOut[["power"]]))
  if(length(targetDelta) == 1 & length(totSampleSizes) > 1)   output$meanPower = matrix(apply(multiThreadOut[["power"]], 2, mean), ncol = 1)
  if(length(targetDelta) > 1)                                 output$meanPower = output$meanPower = apply(multiThreadOut[["power"]], c(2,3), mean)
  rownames(output$meanPower) = totSampleSizes
  colnames(output$meanPower) = targetDelta
  
  # powerArray
  if(length(targetDelta) == 1) output$powerArray = array(data = multiThreadOut[["power"]], dim = c(sims, length(totSampleSizes), length(targetDelta)))
  if(length(targetDelta) > 1 & length(totSampleSizes) == 1) output$powerArray = multiThreadOut[["power"]]
  if(length(targetDelta) > 1) output$powerArray = multiThreadOut[["power"]]
  dimnames(output$powerArray) = list(1:sims, totSampleSizes, targetDelta)  
  
  # deltaArray
  if(length(targetDelta) == 1 & length(totSampleSizes) == 1)  output$deltaArray = list(matrix(multiThreadOut[["delta"]]))
  if(length(targetDelta) == 1 & length(totSampleSizes) > 1)   output$deltaArray = list(multiThreadOut[["delta"]])
  if(length(targetDelta) > 1 & length(totSampleSizes) == 1)   output$deltaArray = lapply(multiThreadOut[["delta"]],as.matrix)
  if(length(targetDelta) > 1 & length(totSampleSizes) > 1)    output$deltaArray = multiThreadOut[["delta"]]
  names(output$deltaArray) = targetDelta
  for(d in 1:length(targetDelta)){
    colnames(output$deltaArray[[d]]) = totSampleSizes
  }
  
  return(output)
}
