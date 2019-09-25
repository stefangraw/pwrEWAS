
#
getAlphBet <- function(myMean, myVar){
    alpha <- myMean^2 * ((1-myMean)/myVar - 1/myMean) 
    beta <- alpha * (1/myMean - 1)
    return( list(alpha = alpha, beta = beta) )
}

#
getMeanVar <- function(myAlpha, myBeta){
    mean <- myAlpha / (myAlpha + myBeta)
    var <- myAlpha * myBeta / ((myAlpha + myBeta)^2 * (myAlpha + myBeta + 1))
    return( list(mean = mean, var = var) )
}

# beta-value to m-value
beta2Mvalue <- function(beta){ # beta to m-value
    return( log2(beta / (1-beta)) )
}

# t-test uses UNequal variance
ttestSlow <- function(g1Beta,g2Beta,rCnt,rTx,paired){
    mvals <- cbind(beta2Mvalue(g1Beta), beta2Mvalue(g2Beta))
    ttest <- apply(mvals,1,  function (x) 
        stats::t.test(x[seq_len(rCnt)], x[(rCnt+1):(rCnt+rTx)], 
            paired = FALSE, 
            var.equal = TRUE))
    temp <- NULL
    temp$pval <- unlist(lapply(ttest, function(x) x$p.value))
    temp$fdr <- stats::p.adjust(temp$pval, method = "fdr")
    return(temp)
}


# t-test uses equal variance
ttestFast <- function(g1Beta,g2Beta,rCnt,rTx){
    mvals <- cbind(beta2Mvalue(g1Beta), beta2Mvalue(g2Beta))
    ttest <- genefilter::rowttests(mvals, 
        fac = factor(c(rep("g1",rCnt),rep("g2",rTx))), 
        tstatOnly = FALSE)  # faster: tstatOnly = T
    temp <- NULL
    temp$pval <- ttest$p.value
    temp$fdr <- stats::p.adjust(temp$pval, method = "fdr")
    return(temp)
}


# Wilcox rank sum test
Wilcox <- function(g1Beta, g2Beta, rCnt,rTx){
    mvals <- cbind(beta2Mvalue(g1Beta), beta2Mvalue(g2Beta))
    WRS <- apply(mvals,1, 
        function (x) stats::wilcox.test(x[seq_len(rCnt)] - 
                x[(rCnt+1):(rCnt+rTx)], correct=TRUE))
    temp <- NULL
    temp$pval <- unlist(lapply(WRS, function(x) x$p.value))
    temp$fdr <- stats::p.adjust(temp$pval, method = "fdr")
    return(temp)
}

# limma
limma <- function(g1Beta, g2Beta, rCnt,rTx){
    mvals <- cbind(beta2Mvalue(g1Beta), beta2Mvalue(g2Beta))
    design <- stats::model.matrix(~ c(rep("g1",rCnt),rep("g2",rTx)))
    limmaFit <- limma::lmFit(mvals, design)
    temp <- NULL
    temp$pval <- eBayes(limmaFit)$p.value[,2]
    temp$fdr <- stats::p.adjust(temp$pval, method = "fdr")
    return(temp)
}

# CPGassoc
CPGassoc <- function(g1Beta, g2Beta, rCnt,rTx){
    assoc <- CpGassoc::cpg.assoc(cbind(g1Beta, g2Beta), 
        c(rep("g1",rCnt),rep("g2",rTx)))
    temp <- NULL
    temp$pval <- assoc$results$P.value
    temp$fdr <- assoc$results$FDR
    return(temp)
}


getTau <- function(targetDmCpGs, targetDelta, methPara, 
    detectionLimit, J, CpGonArray){
    out <- NULL
    tau <- 1
    tauSteps <- 1
    lookForTau <- TRUE
    cnt <- 0
    maxCnt <- 100
    
    while(cnt < maxCnt & lookForTau){
        # print(tau)
        percentile <- NULL
        for(i in seq_len(100)){
            # simulate deltas for J CpG's (number of simulated CpG's later)
            cpgIdx4Tau <- sample(x = seq_len(CpGonArray), 
                size = J, replace = TRUE) # pick J random CpG's to be changed in mean meth
            delta <- truncnorm::rtruncnorm(1, mean = 0, sd = tau, 
                a=0.5 - methPara$mu[cpgIdx4Tau] - 
                    sqrt(0.25-methPara$var[cpgIdx4Tau]), 
                b=0.5 - methPara$mu[cpgIdx4Tau] + 
                    sqrt(0.25-methPara$var[cpgIdx4Tau]))
            # 99.999% percentile 
            percentile[i] <- stats::quantile(abs(delta),0.9999)
        }
        
        # next tau step
        if(mean(percentile) < targetDelta - 0.5*detectionLimit & tau >= 1){
            tau <- tau + 1
        } else if(mean(percentile) < targetDelta - 0.5*detectionLimit & tau < 1){
            tauSteps <- 0.5 * tauSteps
            tau <- tau + tauSteps
        } else if(mean(percentile) > targetDelta + 0.5*detectionLimit){
            tauSteps <- 0.5 * tauSteps
            tau <- tau - tauSteps
        } else {
            lookForTau <- FALSE
        }
        cnt <- cnt + 1
    }
    if(cnt == maxCnt) stop("Max iterations reached")
    truelyDMperc <- mean(abs(delta) > detectionLimit)
    out$tau <- tau
    targetK <- round(1/truelyDMperc * targetDmCpGs)
    out$K <- ifelse(targetK > J, J, targetK)
    return(out)
}




getK <- function(targetDmCpGs, methPara, detectionLimit, J, CpGonArray, tau){
    cpgIdx4Tau <- sample(x = seq_len(CpGonArray), size = J, replace = TRUE) # pick J random CpG's to be changed in mean meth
    delta <- truncnorm::rtruncnorm(1, mean = 0, sd = tau, 
        a=0.5 - methPara$mu[cpgIdx4Tau] - sqrt(0.25-methPara$var[cpgIdx4Tau]), 
        b=0.5 - methPara$mu[cpgIdx4Tau] + sqrt(0.25-methPara$var[cpgIdx4Tau]))
    
    truelyDMperc <- mean(abs(delta) > detectionLimit)
    targetK <- round(1/truelyDMperc * targetDmCpGs)
    K <- ifelse(targetK > J, J, targetK)
    return(K)
}


