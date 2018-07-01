#' @export
myPlotCI3D = function(data#, Nmin, Nmax, Nsteps, deltaSigma
){
  # sampleSizes = seq(Nmin, Nmax, Nsteps)
  sampleSizes = as.numeric(dimnames(data)[[2]])
  deltas = dimnames(data)[[3]]
  
  # df <- data.frame(x = sampleSizes, y = colMeans(data[,,1]))
  df <- data.frame(x = sampleSizes, y = colMeans(matrix(data[,,1])))
  p = ggplot(df, aes(x = x, y = y)) + 
    ggtitle("Mean power curve with 95-percentile interval (2.5% & 97.5%)") + 
    labs(x = "Sample size") + 
    labs(y = "Power") +
    geom_hline(yintercept = 0.8, linetype = 3, size = 1.5)
  dftemp = NULL
  # cols = rainbow(dim(data)[3])
  
  if(length(deltas) == 1){
    scatter = 0
  } else if(length(sampleSizes) > 1 & length(deltas) > 1) {
    scatter = seq(-0.03 * (sampleSizes[2]-sampleSizes[1]), 0.03 * (sampleSizes[2]-sampleSizes[1]), 0.06/(length(deltas)-1) * (sampleSizes[2]-sampleSizes[1]))
  } else scatter = seq(-0.03 , 0.03, 0.06/(length(deltas)-1))
  

    
    
  for(j in 1:dim(data)[3]){
    U = NULL
    L = NULL
    if(class(data[,,j]) == "matrix"){
      dataSlice = data[,,j]
    } else  dataSlice = matrix(data[,,j])
   
    for(i in 1:dim(dataSlice)[2]){
      # error <- qt(0.975,df=length(dataSlice[,i])-1)*sd(dataSlice[,i])/sqrt(length(dataSlice[,i]))
      # L[i] <- mean(dataSlice[,i])-error
      # U[i] <- mean(dataSlice[,i])+error
      L[i] = quantile(dataSlice[,i], 0.025, na.rm = T)
      U[i] = quantile(dataSlice[,i], 0.975, na.rm = T)
    }
    
    dftemp[[j]] <- data.frame(x = sampleSizes, y = colMeans(dataSlice) , 
                              scatter = sampleSizes + scatter[j],
                              L = L, U = U, deltas = rep(as.character(deltas[j]), length(sampleSizes)))
    
    p = p +
      geom_errorbar(data = dftemp[[j]], aes(x = scatter, ymax = U, ymin = L, colour=deltas),
                    width=ifelse(length(sampleSizes)>1,diff(range(sampleSizes))/(length(sampleSizes)*4),0.9), linetype=1, size=1) +
      geom_line(data = dftemp[[j]], aes(x = x, y = y, colour=deltas), size=1.2) +
      geom_point(data = dftemp[[j]], aes(x = x, y = y), size = 1) +
      scale_x_continuous(breaks = sampleSizes) +
      scale_y_continuous(minor_breaks = seq(0 , 1, 0.1), breaks = seq(0, 1, 0.2), limits = c(0,1)) +
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=22)) +
      theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))+ 
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+ 
      theme(legend.text=element_text(size=17),
            legend.title=element_text(size=20)) + 
      scale_colour_discrete(name  = expression(Delta))
    
  }
  print(p)
}
# myPlotCI3D(out$powerArray)

#' @export
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' @export
myDensityPlots = function(data, detectionLimit){
  maxDensY = 0
  maxDensX = 0
  for(d in 1:length(data)){
    dens = density(data[[d]][abs(data[[d]])>detectionLimit]) ## XXXXX need common bandwidth
    maxDensY = max(c(maxDensY, max(dens$y)))
    maxDensX = max(c(maxDensX, max(abs(dens$x))))
  }
  plot(density(data[[1]]), col = "white", ylim = c(0,maxDensY), xlim = c(max(-1,-1.1*maxDensX), min(1,1.1*maxDensX)),
       main = "", xlab = expression(Delta), cex.axis = 1.5, cex.lab = 1.5)
  myLineWd = 2.5
  for(d in 1:length(data)){
    # lines(density(data[[d]][abs(data[[d]])>detectionLimit]), col = gg_color_hue(length(data))[d], lwd = myLineWd)
    lines(density(data[[d]][abs(data[[d]])>detectionLimit], from = detectionLimit), col = gg_color_hue(length(data))[d], lwd = myLineWd)
    lines(density(data[[d]][abs(data[[d]])>detectionLimit], to = -detectionLimit), col = gg_color_hue(length(data))[d], lwd = myLineWd)
  }
  abline(v = c(-detectionLimit, detectionLimit), lty = 3)
}

# myDensityPlots(out$deltaArray, 0.01)
