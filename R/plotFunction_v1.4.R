#' @title Plot funtion to create a power plot
#'
#' @description pwrEWAS_powerPlot create a figure with power (with 95-percentile interval (2.5% & 97.5%)) as a funtion sample size for different effect sizes
#' 
#' @param data "powerArray" attribute within the pwrEWAS object create by pwrEWAS.
#' @param sd FALSE if targetDelta was specified in pwrEWAS, and TRUE if deltaSD was specified in pwrEWAS.
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
#' pwrEWAS_powerPlot(data = outDelta$powerArray, sd = FALSE)
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
#' pwrEWAS_powerPlot(data = outSD$powerArray, sd = TRUE)
pwrEWAS_powerPlot <- function(data, sd = FALSE){
  sampleSizes <- as.numeric(dimnames(data)[[2]])
  deltas <- dimnames(data)[[3]]
  
  df <- data.frame(x = sampleSizes, y = colMeans(matrix(data[,,1])))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) + 
    ggplot2::ggtitle("Mean power curve with 95-percentile interval (2.5% & 97.5%)") + 
    ggplot2::labs(x = "Sample size") + 
    ggplot2::labs(y = "Power") +
    ggplot2::geom_hline(yintercept = 0.8, linetype = 3, size = 1.5)
  dftemp <- NULL
  
  if(length(deltas) == 1){
    scatter <- 0
  } else if(length(sampleSizes) > 1 & length(deltas) > 1) {
    scatter <- seq(-0.03 * (sampleSizes[2]-sampleSizes[1]), 0.03 * (sampleSizes[2]-sampleSizes[1]), 0.06/(length(deltas)-1) * (sampleSizes[2]-sampleSizes[1]))
  } else scatter <- seq(-0.03 , 0.03, 0.06/(length(deltas)-1))
  
  
  
  
  for(j in 1:dim(data)[3]){
    U <- NULL
    L <- NULL
    if(is(data[,,j], "matrix")){
      dataSlice <- data[,,j]
    } else  dataSlice <- matrix(data[,,j])
    
    for(i in 1:dim(dataSlice)[2]){
      L[i] <- quantile(dataSlice[,i], 0.025, na.rm = TRUE)
      U[i] <- quantile(dataSlice[,i], 0.975, na.rm = TRUE)
    }
    
    dftemp[[j]] <- data.frame(x = sampleSizes, y = colMeans(dataSlice) , 
                              scatter = sampleSizes + scatter[j],
                              L = L, U = U, deltas = rep(as.character(deltas[j]), length(sampleSizes)))
    
    p <- p +
      ggplot2::geom_errorbar(data = dftemp[[j]], ggplot2::aes(x = scatter, ymax = U, ymin = L, colour=deltas),
                             width=ifelse(length(sampleSizes)>1,diff(range(sampleSizes))/(length(sampleSizes)*4),0.9), linetype=1, size=1) +
      ggplot2::geom_line(data = dftemp[[j]], ggplot2::aes(x = x, y = y, colour=deltas), size=1.2) +
      ggplot2::geom_point(data = dftemp[[j]], ggplot2::aes(x = x, y = y), size = 1) +
      ggplot2::scale_x_continuous(breaks = sampleSizes) +
      ggplot2::scale_y_continuous(minor_breaks = seq(0 , 1, 0.1), breaks = seq(0, 1, 0.2), limits = c(0,1)) +
      ggplot2::theme(axis.text=ggplot2::element_text(size=18),
                     axis.title=ggplot2::element_text(size=22)) +
      ggplot2::theme(axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0)))+ 
      ggplot2::theme(axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0)))+ 
      ggplot2::theme(legend.text=ggplot2::element_text(size=17),
                     legend.title=ggplot2::element_text(size=20))
    
    if(sd){
      p <- p + ggplot2::scale_colour_discrete(name  = expression(paste("sd(",Delta[beta],")",sep = "")))
    } else p <- p + ggplot2::scale_colour_discrete(name  = expression(Delta[beta]))
    
    
  }
  print(p)
}


gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' @title Density plot for simulated differences in mean methylation
#'
#' @description pwrEWAS_deltaDensity create a density plot of the simulated differences in mean methylation for different effect sizes
#' 
#' @param data "deltaArray" attribute within the pwrEWAS object create by pwrEWAS
#' @param detectionLimit Detection limit specified in pwrEWAS.
#' @param sd FALSE if targetDelta was specified in pwrEWAS, and TRUE if deltaSD was specified in pwrEWAS.
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
#' pwrEWAS_deltaDensity(data = outDelta$deltaArray, detectionLimit = 0.01, sd = FALSE)
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
#' pwrEWAS_deltaDensity(data = outSD$deltaArray, detectionLimit = 0.01, sd = TRUE)
pwrEWAS_deltaDensity <- function(data, detectionLimit = 0.01, sd = FALSE){
  maxDensY <- 0
  maxDensX <- 0
  for(d in 1:length(data)){
    dens <- density(data[[d]][abs(data[[d]])>detectionLimit]) ## XXXXX need common bandwidth
    maxDensY <- max(c(maxDensY, max(dens$y)))
    maxDensX <- max(c(maxDensX, max(abs(dens$x))))
  }
  plot(density(data[[1]]), col = "white", ylim = c(0,maxDensY), xlim = c(max(-1,-1.1*maxDensX), min(1,1.1*maxDensX)),
       main = "", xlab = expression(Delta[beta]), cex.axis = 1.5, cex.lab = 1.5)
  myLineWd <- 2.5
  for(d in 1:length(data)){
    lines(density(data[[d]][abs(data[[d]])>detectionLimit], from = detectionLimit), col = gg_color_hue(length(data))[d], lwd = myLineWd)
    lines(density(data[[d]][abs(data[[d]])>detectionLimit], to = -detectionLimit), col = gg_color_hue(length(data))[d], lwd = myLineWd)
  }
  abline(v = c(-detectionLimit, detectionLimit), lty = 3)
  if(sd){
    legend("topright", names(data), col = gg_color_hue(length(data)), lty = 1, lwd = myLineWd, title = expression(paste("sd(",Delta[beta],")",sep = "")), cex = 1.5)
  } else legend("topright", names(data), col = gg_color_hue(length(data)), lty = 1, lwd = myLineWd, title = expression(Delta[beta]), cex = 1.5)
  
}


