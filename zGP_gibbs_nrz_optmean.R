zGP_gibbs_nrz_optmean = function(xd, yRL, Ngibbs, locs, Nzfit){
  locx = locs[1]
  locy = locs[2]
  y = yRL
  
  N = length(y)
  Nall = N
  Ndim = dim(xd)[2]
  
  inds = which(y <= 0)
  xz = xd[inds,]
  indsz = inds
  Nz = length(indsz)
  yz = matrix(y[indsz,])
  
  indsp = which(y > 0)
  Np = length(indsp)
  xp = xd[indsp,]
  yp = matrix(y[indsp,])
  
  ypp = matrix(y[1:(Np+Nzfit),])
  xpp = xd[1:(Np+Nzfit),]
  xall = rbind(xz, xp)
  yall = matrix(rbind(yz, yp))
  
  sigstarsave = matrix(0, nrow = length(indsz), ncol = Ngibbs)
  ysmp = yall
  ykeep = ysmp
  yorg = ykeep # same as y, but order negatives to positives -- goes with xall
  ysmpsave = matrix(0, nrow = Nall, ncol = Ngibbs)
  rhosqsave = matrix(0, nrow = Ngibbs, ncol = dim(xall)[2])
  
  mu = function(x) cbind(matrix(1, nrow = dim(x)[1], ncol = 1), x)%*%B # there is no B defined before?
  
  mat52 = function(d) (1+sqrt(5)*d+(5/3)*d^2)*exp(-sqrt(5)*d)
  
  ysmp = yall
  if (sum(locs) == 0){
    meantrend = 'linear_mean'
  } else {
    meantrend = 'physical_mean'
  }
  switch(meantrend, # you can define your own trend here if you'd like
         'zero_mean' = {
           # B = mean(ysmp)
           # mu = function(x) matrix(1, nrow = dim(x)[1], ncol = 1)%*%B
           ppgasp_options$zero_mean = "Yes"
           ppgasp_options$nugget.est = F
           model = ppgasp(design = xpp, response = ypp - mean(ypp),
                          zero.mean = ppgasp_options$zero.mean,
                          nugget.est = ppgasp_options$nugget.est)
         },
         'physical_mean' = {
           # H = function(xd) cbind(matrix(1, nrow = dim(xd)[1], ncol = 1), matrix(xd[,1]), matrix(sqrt((xd[,3]-locx)^2 + (xd[,4]-locy)^2)))
           # 1- vol, 2- BF, 3- E, 4- N -- specific for Aluto
           # Probably should be
           # H = function(xd) cbind(matrix(1, nrow = dim(xd)[1], ncol = 1), matrix(xd[,1]), matrix(xd[,2]), matrix(sqrt((xd[,3]-locx)^2 + (xd[,4]-locy)^2)))
           # B = solve(t(H(xall))%*% H(xall)) %*% t(H(xall)) %*% yall
           # mu = function(x) H(x)%*%B
           ppgasp_options$trend = cbind(matrix(1, nrow = Np+Nzfit, ncol = 1), xpp[,1], xpp[,2], sqrt((xpp[,4]-locx)^2+(xpp[,5]-locy)^2))
           ppgasp_options$zero.mean = "No"
           ppgasp_options$nugget.est = F
           model = ppgasp(design = xpp, response = ypp,
                          trend = ppgasp_options$trend,
                          zero.mean = ppgasp_options$zero.mean,
                          nugget.est = ppgasp_options$nugget.est)
         },
         'linear_mean' = {
           ppgasp_options$trend = cbind(matrix(1, nrow = Np+Nzfit, ncol = 1), xpp)
           ppgasp_options$zero.mean = "No"
           ppgasp_options$nugget.est = F
           model = ppgasp(design = xpp, response = ypp,
                          trend = ppgasp_options$trend,
                          zero.mean = ppgasp_options$zero.mean,
                          nugget.est = ppgasp_options$nugget.est)
         })
  
  thetas1 = 1/model@beta_hat
  rhosq = t(thetas1^2)
  
  ykeep = ysmp
  ysmpgr = ysmp
  yorg = ykeep
  count = 0
  
  Ball = corr_matrix(xall, xall, rhosq)
  Ball = Ball + (1e-6)*diag(dim(Ball)[1]) #adding the nugget 1e-3 too big, 1e-12 too small
  invBall = solve(Ball)
  sigs = as.numeric((1/Nall)*t(ysmp)%*%(solve(Ball)%*%ysmp))
  
  sigssave = c()
  
  for (ii in 1:Ngibbs){
    if (ii%%500 == 0){
      my_print = sprintf("ii = %d", ii); print(my_print)
    } # below sample range pars every 50th step
    if (ii == 1 | ii%%50 == 0){
      if (ii == 1){
        ypp = matrix(y[1:(Np+Nzfit)])
      } else {
        ypp = matrix(ysmpgr[1:(Np+Nzfit)])
      }
      model = ppgasp(design = xpp, response = ypp,
                     trend = ppgasp_options$trend,
                     zero.mean = ppgasp_options$zero.mean,
                     nugget.est = ppgasp_options$nugget.est) # ypps change so need to find new range pars and B
      B = model@theta_hat
      thetas1 = 1/model@beta_hat
      rhosq = t(thetas1^2)
      switch(meantrend, # you can define your own trend here if you'd like
             'zero_mean' = {
               B = mean(ysmp)
               mu = function(x) matrix(1, nrow = dim(x)[1], ncol = 1)%*%B
             },
             'physical_mean' = {
               H = function(xpp) cbind(matrix(1, nrow = dim(xpp)[1], ncol = 1), xpp[,1], xpp[,2], sqrt((xpp[,4]-locx)^2+(xpp[,5]-locy)^2))
               # 1- vol, 2- BF, 3- E, 4- N -- specific for Aluto
               # Probably should be
               # H = function(xd) cbind(matrix(1, nrow = dim(xd)[1], ncol = 1), matrix(xd[,1]), matrix(xd[,2]), matrix(sqrt((xd[,3]-locx)^2 + (xd[,4]-locy)^2)))
               # B = solve(t(H(xall))%*% H(xall)) %*% t(H(xall)) %*% yall
               mu = function(x) H(x)%*%B
             },
             'linear_mean' = {
               H = function(x) cbind(matrix(1, nrow = dim(x)[1], ncol = 1), x)
               mu = function(x) H(x)%*%B
             })
      
      Ball = corr_matrix(xpp, xpp, rhosq)
      Ball = Ball + (1e-6)*diag(dim(Ball)[1])
      sigs = as.numeric((1/(Np+Nzfit))*t(ypp)%*%(solve(Ball)%*%ypp))
      
      Dsave = matrix(0, nrow = Nall - 1, ncol = Nz)
      Ballinvsavematrix = matrix(0, nrow = Nall - 1, ncol = Nall - 1)
      Ballinvsave = array(rep(Ballinvsavematrix, Nz), dim = c(Nall - 1, Nall - 1, Nz))
      
      Ball = corr_matrix(xall, xall, rhosq)
      
      Ball = Ball + (1e-6)*diag(dim(Ball)[1])
      invBall = solve(Ball)
      
      Ballinv = solve(Ball)
      Ballinvkeep = Ballinv
      
      D = matrix(1, nrow = Nall - 1, ncol = 1)
      for (kk in 1:Nz){
        Balltemp = Ball[,-kk]
        Balltemp = Balltemp[-kk,]
        xalltemp = xall[-kk,]
        for (jj in 1:(Nall-1)){
          d = sqrt(((xall[kk,] - xalltemp[jj,])^2)/rhosq)
          D[jj,1] = prod(mat52(d))
        }
        Dsave[,kk] = D
        Ballinv = solve(Balltemp)
        Ballinvsave[,,kk] = Ballinv
      }
    }
    ###########################################################################
    # We'll think about one Gibb's sample as a pass across all 0's where we draw
    # one sample at a time conditioned on all of the other samples. This process
    # ends up will correlated (pass to pass) samples so we will only keep every
    # 5th sample below.
    ###########################################################################
    
    for (jj in 1:Nz){
      rr = sample(1:Nz, 1)
      ysmp = ykeep
      ys = matrix(ysmp[-rr,])
      xs = xall[-rr,]
      
      D = matrix(Dsave[,rr])
      Binv = Ballinvsave[,,rr]
      sigsqstar = sigs*(1 - as.numeric(t(D)%*%Binv%*%D))
      sigstarsave[rr,ii] = sigsqstar
      mustar = as.numeric(mu(xall[rr,, drop = F])) + as.numeric(t(D)%*%(Binv%*%(ys - mu(xs))))
      
      # inverse CDF version of sampling
      yatzero = 1/2*(1 + erf(-mustar/(sqrt(2)*sigsqstar)))
      Y = yatzero*runif(1)
      ypsamp = sqrt(2)*sigsqstar*erfinv(2*Y - 1) + mustar
      if (is.infinite(ypsamp)){
        ypsamp = -runif(1)/100
      }
      yg = ypsamp
      ykeep[rr] = yg
    }
    sigssave[ii] = sigs
    ysmp = ykeep
    ysmpgr[indsz] = matrix(ysmp[1:Nz])
    ysmpgr[indsp] = matrix(ysmp[-(1:Nz)])
    ysmpsave[,ii] = ykeep
    rhosqsave[ii,] = rhosq
  }
  ygsmps = ysmpsave
  rhogsmps = rhosqsave
  
  xdr = xd
  xdr[indsz,] = xall[1:Nz,]
  xdr[indsp,] = xall[-(1:Nz),]
  ygr = ygsmps
  ygr[indsz,] = ygsmps[1:Nz,]
  ygr[indsp,] = ygsmps[-(1:Nz),]
  
  output = list(ygr, sqrt(rhogsmps), sigssave)
}
