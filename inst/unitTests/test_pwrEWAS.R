test_pwrEWAS <- function(){
    checkEqualsNumeric(pwrEWAS(minTotSampleSize = 10, 
                               maxTotSampleSize = 10, 
                               SampleSizeSteps = 1, 
                               NcntPer = 0.5, 
                               targetDelta = 0.8, 
                               J = 10000,
                               targetDmCpGs = 100, 
                               tissueType = "Adult (PBMC)", 
                               detectionLimit = 0.01, 
                               DMmethod = "limma", 
                               FDRcritVal = 0.05, 
                               core = 1, 
                               sims = 50)$meanPower,
                       0.6711375, 
                       tolerance = 0.05)
}
