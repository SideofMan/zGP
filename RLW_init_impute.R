RLW_init_impute = function(xd,y){
  N = length(y)
  inds = which(ystart == 0)
  xz = xd[inds,]
  indsz = inds
  Nz = length(indsz)
  
  inds = which(y > 0)
  xp = xd[inds,]
  yp = matrix(y[inds])
  indsp = inds
  Np = length(indsp)
  ysave = y
  
  H = function(xd) matrix(0, nrow = dim(xp)[1], ncol = 1)
  
  model = ppgasp(design = xp, response = yp - mean(yp))
  thetas = 1/model@beta_hat # betas are inverse of range parameters
  
  rhosq = t(thetas^2)
  
  mat52 = function(d) (1+sqrt(5)*d+(5/3)*d^2)*exp(-sqrt(5)*d)
  
  ##############################################
  # This is where we start the negative samples
  ##############################################
  
  yimputesave = matrix(0, nrow = N, ncol = 50)
  for (ii in 1:50){
    print(ii)
    y = ysave
    if (ii == 1){
      B = corr_matrix(xp, xp, rhosq)
      B = B + (1e-6)*diag(dim(B)[1]) # adding the nugget (~0) to the
                                      # diagonal of the matrix
      Btemp = corr_matrix(xd, xd, rhosq)
      Btemp[-(1:Np),] = 0
      Btemp[,-(1:Np)] = 0
      
      Btemp = Btemp + (1e-6)*diag(dim(Btemp)[1]) # adding the nugget (~0) to the
                                                  # diagonal of the matrix
      saveB = B
      saveBtemp=Btemp
      invB=solve(B)
      saveBinv=invB
      invBtemp=solve(Btemp)
      saveBtempinv=invBtemp
    } else {
      B = saveB
      Btemp = saveBtemp
      invB = saveBinv
      invBtemp = saveBtempinv
    }

    # define temp variables for anonymous functions so they are not
    # change later on when the anonymous function is called
    yp_temp1 = yp
    #yc_temp = yc
    
    xp_temp1 = xp
    rhosq_temp1 = rhosq
    Np_temp1 = Np
    invB_temp1 = invB
    
    r = function(z) mat52(sqrt((matrix(rep(z, times = dim(xp_temp1)[1]), nrow = dim(xp_temp1)[1], byrow = T)-xp_temp1)^2%*%t(1/rhosq_temp1)))
    ypred = function(z) t(r(z))%*%(invB_temp1%*%yp_temp1)
    sigs = as.numeric((1/Np)*t(yp)%*%(invB%*%yp))
    
    Cz = corr_matrix(xz, xz, rhosq) # correlation between censored points (C in doc)
    Czo = corr_matrix(xp, xz, rhosq) # correlation between censored points and positive points (D in doc)
    
    # better to define ypp here
    ypp = matrix(0, nrow = 1, ncol = Np)
    for (k in 1:Np){
      ypp[1,k] = ypred(xp[k,])
    }
    mustar = Czo%*%(invB%*%(t(ypp)))
    
    #################################
    # Initial smps to start RLW alg.
    #################################
    
    SigSqstar = sigs*(Cz - Czo%*%(invB%*%t(Czo)))
    Sigma = (SigSqstar+t(SigSqstar))/2
    
    smpsize = 1e4 # defines the sample size from the multivariate Gaussian
                  # (with #variables equal to Nc)
    
    smps = mvrnorm(n = smpsize, mu = mustar, Sigma = Sigma)
    
    isneg = matrix(0, nrow = smpsize, ncol = Nz)
    
    ###########################################################
    # Batch algorithm for initializing negative-4-zeros sample
    ###########################################################
    
    for (k in 1:Nz){
      LL = 0
      indsL = which(smps[,k] < LL)
      isneg[indsL, k] = 1
    }
    inds = order(rowSums(isneg), decreasing = T)
    
    yc = smps[inds[1],]
    inds = which(yc < 0)
    y[indsz[inds]] = yc[inds]
    
    inds = which(y == 0)
    indspp = setdiff(1:N, inds)
    Npp = length(indspp)
    Ncc = length(inds)
    
    indscc = inds
    xpp = xd[indspp,]
    xcc = xd[indscc,]
    ypp = matrix(y[indspp])
    ycc = matrix(y[indscc])
    count = 0
    while (Npp < N){
      count = count + 1
      B = corr_matrix(xpp, xpp, rhosq)
      
      B = B + (1e-6)*diag(dim(B)[1]) # adding the nugget (~0) to the
                                      # diagonal of the matrix
      
      # define temp variables for anonymous functions so they are not
      # change later on when the anonymous function is called
      ypp_temp2 = matrix(ypp)
      ycc_temp2 = matrix(ycc)
      
      xpp_temp2 = xpp
      rhosq_temp2 = rhosq
      Np_temp2 = Np
      B_temp2 = B
      
      solB = solve(B)
      
      r = function(z) mat52(sqrt((matrix(rep(z, times = dim(xpp_temp2)[1]), nrow = dim(xpp_temp2)[1], byrow = T)-xpp_temp2)^2%*%t(1/rhosq_temp2)))
      ypred = function(z) t(r(z))%*%(solB%*%ypp_temp2)
      sigs = as.numeric((1/Npp)*t(ypp)%*%(solB%*%ypp))
      
      if (length(indscc) == 1){
        xcc = matrix(xd[indscc,], nrow = 1)
      }
      Cz = corr_matrix(xcc, xcc, rhosq)
      Czo = corr_matrix(xpp, xcc, rhosq)
      
      # better to define ypp here
      ypp = matrix(0, nrow = 1, ncol = Npp) # this is defined the transpose of what it should be?
      for (k in 1:Npp){
        ypp[1,k] = ypred(xpp[k,])
      }
      
      mustar = Czo%*%(solB%*%t(ypp))
      SigSqstar = sigs*(Cz - Czo%*%(solB%*%t(Czo)))
      Sigma = (SigSqstar+t(SigSqstar))/2
      
      smpsize = 1e5
      smps = mvrnorm(n = smpsize, mu = mustar, Sigma = Sigma)
      isneg = matrix(0, nrow = smpsize, ncol = Ncc)
      for (k in 1:Ncc){
        LL = 0
        indsL = which(smps[,k] < LL)
        isneg[indsL, k] = 1
      }
      
      inds = order(rowSums(isneg), decreasing = T)
      
      yctemp = smps[inds[1],]
      
      inds = which(yctemp < 0)
      
      y[indscc[inds]] = yctemp[inds]
      
      if (count > 50){
        y[indscc] = -runif(length(indscc))/100
      }
      
      inds = which(y == 0)
      indspp = setdiff(1:N, inds)
      Npp = length(indspp)
      Ncc = length(inds)
      indscc = inds
      xpp = xd[indspp,]
      xcc = xd[indscc,] 
      ypp = y[indspp]
      ycc = y[indscc]
    }
    yimputestart = y
    yimputesave[,ii] = yimputestart
  }
  output = list(yimputesave, sigs)
  return(output)
}