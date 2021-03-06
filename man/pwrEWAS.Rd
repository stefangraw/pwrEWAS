% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pwrEWAS_v1.8.R
\name{pwrEWAS}
\alias{pwrEWAS}
\title{pwrEWAS - A computationally efficient tool for comprehensive power estimation in EWAS}
\usage{
pwrEWAS(minTotSampleSize, maxTotSampleSize, SampleSizeSteps, NcntPer,
    targetDelta = NULL, deltaSD = NULL, J = 1e+05, targetDmCpGs,
    tissueType = c("Adult (PBMC)", "Saliva", "Sperm", "Lymphoma",
    "Placenta", "Liver", "Colon", "Blood adult", "Blood 5 year olds",
    "Blood newborns", "Cord-blood (whole blood)", "Cord-blood (PBMC)"),
    detectionLimit = 0.01, DMmethod = c("limma", "t-test (unequal var)",
    "t-test (equal var)", "Wilcox rank sum", "CPGassoc"),
    FDRcritVal = 0.05, core = 1, sims = 50)
}
\arguments{
\item{minTotSampleSize}{Minimum total sample size.}

\item{maxTotSampleSize}{Maximum total sample size.}

\item{SampleSizeSteps}{Sample size increments.}

\item{NcntPer}{Percentage sample group 1 (control group) (NcntPer = 0.5 indicates a balanced design).}

\item{targetDelta}{Target maximum difference in mean DNAm. (Either 'targetDelta' or 'deltaSD' should be specified)}

\item{deltaSD}{Standard deviation of simulated differences. (Either 'targetDelta' or 'deltaSD' should be specified)}

\item{J}{Number of CpGs tested/simulated (default: 100000).}

\item{targetDmCpGs}{Target number of DM CpGs.}

\item{tissueType}{Select a tissue type from the list of most commonly used tissue types: "Adult (PBMC)" (default), "Saliva", "Sperm", "Lymphoma", "Placenta", "Liver", "Colon", "Blood adult", "Blood 5 year olds", "Blood newborns", "Cord-blood (whole blood)" or "Cord-blood (PBMC)".}

\item{detectionLimit}{Smallest detectable difference in DNAm (default: 0.01).}

\item{DMmethod}{Method of Differential Methylation analysis: "limma" (default), "t-test (unequal var)", "t-test (equal var)", "Wilcox rank sum", "CPGassoc".}

\item{FDRcritVal}{FDRcritVal (default: 0.05).}

\item{core}{Number of threads for multi-threading (default: 1).}

\item{sims}{Number of simulated data sets (default: 50).}
}
\value{
pwrEWAS will return an object with the following four attributes: meanPower, powerArray, deltaArray, and metric, where metric contains marTypeI, classicalPower, FDR, and FDC
}
\description{
pwrEWAS is a computationally efficient tool to estimate power in EWAS as a function of sample and effect size for two-group comparisons of DNAm (e.g., case vs control, exposed vs non-exposed, etc.). Detailed description of in-/outputs, instructions and an example, as well as interpretations of the example results are provided in the vignette: vignette("pwrEWAS") 
}
\examples{
outDelta <- pwrEWAS(minTotSampleSize = 10,
    maxTotSampleSize = 20,
    SampleSizeSteps = 10,
    NcntPer = 0.5,
    targetDelta = c(0.2, 0.5),
    J = 1000,
    targetDmCpGs = 10,
    tissueType = "Adult (PBMC)",
    detectionLimit = 0.01,
    DMmethod = "limma",
    FDRcritVal = 0.05,
    core = 2,
    sims = 30)
    
outSD <- pwrEWAS(minTotSampleSize = 10,
    maxTotSampleSize = 20,
    SampleSizeSteps = 10,
    NcntPer = 0.5,
    deltaSD = c(0.02, 0.03),
    J = 1000,
    targetDmCpGs = 10,
    tissueType = "Adult (PBMC)",
    detectionLimit = 0.01,
    DMmethod = "limma",
    FDRcritVal = 0.05,
    core = 2,
    sims = 30)
}
\keyword{DNAm}
\keyword{microarray}
\keyword{power}
