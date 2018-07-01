
#' @export
getAlphBet = function(myMean, myVar){
  alpha = myMean^2 * ((1-myMean)/myVar - 1/myMean) 
  beta = alpha * (1/myMean - 1)
  return( list(alpha = alpha, beta = beta) )
}

#' @export
getMeanVar = function(myAlpha, myBeta){
  mean = myAlpha / (myAlpha + myBeta)
  var = myAlpha * myBeta / ((myAlpha + myBeta)^2 * (myAlpha + myBeta + 1))
  return( list(mean = mean, var = var) )
}

# beta-value to m-value
#' @export
beta2Mvalue = function(beta){ # beta to m-value
  # # g1Mval = log2(g1Beta/(1-g1Beta))
  # # g2Mval = log2(g2Beta/(1-g2Beta))
  # # library(car)
  # g1Mval = logit(g1Beta)
  # g2Mval = logit(g2Beta)
  # # library(lumi)
  # # g1Mval = beta2m(g1Beta)
  # # g2Mval = beta2m(g2Beta)
  # return(suppressWarnings(logit(g1Beta)))
  return( log2(beta / (1-beta)) )
}

# t-test uses UNequal variance
#' @export
ttestSlow = function(g1Beta,g2Beta,rCnt,rTx,paired){
  mvals = cbind(beta2Mvalue(g1Beta), beta2Mvalue(g2Beta))
  ttest = apply(mvals,1,  function (x) t.test(x[1:rCnt], x[(rCnt+1):(rCnt+rTx)], paired = pairedSamples, var.equal = T))
  temp = NULL
  temp$pval = unlist(lapply(ttest, function(x) x$p.value))
  temp$fdr = p.adjust(temp$pval, method = "fdr")
  # ttest_u_alpha = sum(ttest_u_pval<0.05)/J
  return(temp)
}


# t-test uses equal variance
#' @export
ttestFast = function(g1Beta,g2Beta,rCnt,rTx){
  mvals = cbind(beta2Mvalue(g1Beta), beta2Mvalue(g2Beta))
  ttest = rowttests(mvals, fac = factor(c(rep("g1",rCnt),rep("g2",rTx))), tstatOnly = F)  # faster: tstatOnly = T
  temp = NULL
  temp$pval = ttest$p.value
  temp$fdr = p.adjust(temp$pval, method = "fdr")
  # ttest_u_alpha = sum(ttest_u_pval<0.05)/J
  return(temp)
}


# Wilcox rank sum test
#' @export
Wilcox = function(g1Beta, g2Beta, rCnt,rTx){
  mvals = cbind(beta2Mvalue(g1Beta), beta2Mvalue(g2Beta))
  WRS = apply(mvals,1,  function (x) wilcox.test(x[1:rCnt] - x[(rCnt+1):(rCnt+rTx)], correct=T))
  temp = NULL
  temp$pval = unlist(lapply(WRS, function(x) x$p.value))
  # temp$alpha = sum(WRS_pval<0.05)/J
  temp$fdr = p.adjust(temp$pval, method = "fdr")
  return(temp)
}

# limma
#' @export
limma = function(g1Beta, g2Beta, rCnt,rTx){
  mvals = cbind(beta2Mvalue(g1Beta), beta2Mvalue(g2Beta))
  design <- model.matrix(~ c(rep("g1",rCnt),rep("g2",rTx)))
  limmaFit <- lmFit(mvals, design)
  temp = NULL
  temp$pval <- eBayes(limmaFit)$p.value[,2]
  # limmaFit_alpha <- sum(limmaFit_pval<0.05)/J
  temp$fdr = p.adjust(temp$pval, method = "fdr")
  return(temp)
}

# CPGassoc
#' @export
CPGassoc = function(g1Beta, g2Beta, rCnt,rTx){
  assoc = cpg.assoc(cbind(g1Beta, g2Beta), c(rep("g1",rCnt),rep("g2",rTx)))
  temp = NULL
  temp$pval = assoc$results$P.value
  # temp$alpha = sum(assoc_pval<0.05)/J
  temp$fdr = assoc$results$FDR
  return(temp)
}


#' @export
getTau = function(targetDmCpGs, targetDelta, methPara, detectionLimit, J, CpGonArray){
  library("truncnorm")
  out = NULL
  tau = 1
  tauSteps = 1
  lookForTau = TRUE
  cnt = 0
  maxCnt = 100
  
  while(cnt < maxCnt & lookForTau){
    # print(tau)
    percentile = NULL
    for(i in 1:100){
      # simulate deltas for J CpG's (number of simulated CpG's later)
      cpgIdx4Tau = sample(x = 1:CpGonArray, size = J, replace = T) # pick K random CpG's to be changed in mean meth
      delta = rtruncnorm(1, mean = 0, sd = tau, 
                         a=0.5 - methPara$mu[cpgIdx4Tau] - sqrt(0.25-methPara$var[cpgIdx4Tau]), 
                         b=0.5 - methPara$mu[cpgIdx4Tau] + sqrt(0.25-methPara$var[cpgIdx4Tau]))
      # 99% percentile 
      percentile[i] = quantile(abs(delta),0.99)
    }
    # print(mean(percentile))
    
    # next tau step
    if(mean(percentile) < targetDelta - 0.5*detectionLimit & tau >= 1){
      tau = tau + 1
    } else if(mean(percentile) < targetDelta - 0.5*detectionLimit & tau < 1){
      tauSteps = 0.5 * tauSteps
      tau = tau + tauSteps
    } else if(mean(percentile) > targetDelta + 0.5*detectionLimit){
      tauSteps = 0.5 * tauSteps
      tau = tau - tauSteps
    } else {
      lookForTau = FALSE
    }
    cnt = cnt + 1
  }
  if(cnt == maxCnt) stop("Max iterations reached")
  truelyDMperc = mean(abs(delta) > detectionLimit)
  out$tau = tau
  out$K = round(1/truelyDMperc * targetDmCpGs)
  return(out)
}

#' @export
loadDataset = function(tissueType){
  if(tissueType == "Saliva") load("dataSets/GSE92767.Rdata") else 
  if(tissueType == "Lymphoma") load("dataSets/GSE42372.Rdata") else 
  if(tissueType == "Placenta") load("dataSets/GSE62733.Rdata") else 
  if(tissueType == "Liver") load("dataSets/GSE61258.Rdata") else 
  if(tissueType == "Colon") load("dataSets/GSE77718.Rdata") else 
  if(tissueType == "Peripheral Leukocytes") load("dataSets/GSE42861.Rdata") else
  if(tissueType == "Blood 5 year olds") load("dataSets/GSE83334.Rdata") else 
  if(tissueType == "Blood newborns") load("dataSets/GSE82273.Rdata") else
  if(tissueType == "Cord-blood (whole blood)") load("dataSets/GSE69176.Rdata") else
  if(tissueType == "Cord-blood (PBMC)") load("dataSets/GSE110128.Rdata") else
  if(tissueType == "Adult (PBMC)") load("dataSets/GSE67170.Rdata") else
  stop("Tissue type not found")

  return(methPara)
}



