probs_zeros = function(Nz, Np, xdp, yp, xdn, yRL, indsp, indsz){
  mat52 = function(d) (1+sqrt(5)*d+(5/3)*d^2)*exp(-sqrt(5)*d)
  Np = length(indsp)
  B = matrix(0, nrow = Nz, ncol = Nz) # correlation matrix of xp points
  Btemp = matrix(0, nrow = Nz, ncol = Nz)
  
  # new version of code where ppgasp() does linear trend for you
  
  ppgasp_options = list("trend" = cbind(matrix(1, nrow = Np, ncol = 1), xdp),
                        "zero.mean" = "No",
                        "nugget.est" = F,
                        "testing.trend" = cbind(matrix(1, nrow = Nz, ncol = 1), xdn))
  
  modelp = ppgasp(xdp, yp,
                  trend = ppgasp_options$trend,
                  zero.mean = ppgasp_options$zero.mean,
                  nugget.est = ppgasp_options$nugget.est)
  
  pred_model = predict.ppgasp(modelp, xdn,
                              testing_trend = ppgasp_options$testing.trend)
  
  pmean = pred_model$mean
  psd = pred_model$sd
  
  ####################################################
  # plot of figure 1
  
  plot(1:length(pmean), pmean,
       xlim = c(0,length(pmean)), ylim = c(floor(min(pmean)-2),ceiling(max(pmean)+2)), xaxs = 'i',
       axes = F,
       xlab = "", ylab = "",
       pch = "*", col = "blue")
  axis(1, at = seq(0,length(pmean), by = 10), tck = 0.01)
  axis(2, at = seq(floor(min(pmean)-2),ceiling(max(pmean)+2), by = 2), tck = 0.01)
  box()
  
  erf = function(x) 2*pnorm(sqrt(2)*x)-1
  erfinv = function (x) qnorm((1 + x)/2)/sqrt(2)
  
  Pn = c()
  for (k in 1:length(pmean)){
    lines(c(k,k), c(pmean[k]-2*psd[k], pmean[k]+2*psd[k]), col = "blue")
    Pn[k] = 1/2*(1 + erf(-pmean[k]/(psd[k]*sqrt(2))))
  }
  lines(c(0,Nz), c(0,0), lwd = 3, col = "blue")
  ###################################################
  
  for (j in 1:Nz){
    for (k in 1:Np){
      distnp[j,k] = sqrt((xdn[j,1] - xdp[k,1])^2 + (xdn[j,2] - xdp[k,2])^2)
    }
  }
  
  mindist = matrix(apply(t(distnp), MARGIN = 2, min))
  indsmp = order(mindist)
  indspn = order(Pn)
  Ntemp = round(0.66*length(indspn)) # 0.66 is chosen arbitrarily (about 2/3 of each list)
  indsadd = sort(intersect(indspn[1:Ntemp], indsmp[1:Ntemp]))
  
  Ninclude = length(indsadd)
  indsrest = setdiff(1:Nz, indsadd)
  xzs = rbind(xdn[indsadd,], xdn[indsrest,])
  xall = rbind(xdp, xzs)
  yall = rbind(matrix(yRL[indsp]), matrix(yRL[indsz[indsadd]]), matrix(yRL[indsz[indsrest]]))
  xp_zp = rbind(xdp, xdn[indsadd,])
  yp_zp = rbind(matrix(yRL[indsp]), matrix(yRL[indsz[indsadd]]))
  
  ppgasp_options$trend = cbind(matrix(1, nrow = Np+Ninclude, ncol = 1), xp_zp)
  
  model = ppgasp(xp_zp, yp_zp,
                 trend = ppgasp_options$trend,
                 zero.mean = ppgasp_options$zero.mean,
                 nugget.est = ppgasp_options$nugget.est)
  
  #################################################################
  # plot of figure 2
  plot(mindist, Pn, pch = "*", col = "blue", xaxs = 'i', yaxs = 'i', axes = F,
       xlim = c(0,max(mindist)), ylim = c(min(Pn),1))
  x_ticks = seq(0,max(mindist), by = 0.1)
  y_ticks = seq(0.05*floor(1/0.05*min(Pn)),1, by = 0.05)
  axis(1, at = x_ticks, labels = x_ticks, tck = 0.01)
  axis(2, at = y_ticks, labels = y_ticks, tck = 0.01)
  box()
  points(mindist[indsadd,], Pn[indsadd], pch = "*", col = "red")
  #################################################################
  output = list(xall, yall, Ninclude, ppgasp_options, erf, erfinv)
}