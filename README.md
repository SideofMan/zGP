## R code for zero censored Gaussian process emulation

Function descriptions:

zGP:

The main script that runs the whole algorithm.  
It takes in design coordinates and outputs as xd and y, and runs the zero-censored Gaussian process

RLW_init_impute (Algorithm 2):

This function initializes negative imputed samples for the zero censored points.

zGP_gibbs_nrz_optmean (Algorithm 1):

This function performs zGP substitution sampling to generate a sequence of negative imputed values for the zero-censored points in the design,
whose distribution is converging to the desired Gaussian Process on the design conditional that the GP's outputs are equal to the design's outputs where the outputs are positive, and the GP's outputs are negative where the design's outputs are zero.

corr_matrix:

This function is a helper function that creates a correlation matrix between 2 sets of points (eqn. 2.1).

probs_zeros:

This function calculates the probability of initial negative imputed values of being zeros and determines which zero-censored points in the design should be used to fit the GP along with the design points whose outputs are positive (section 3.3.2).

Keys for some variable names:

zGP

xd: design coordinates
ytrue: true outputs from toy function  
Nd: # of points  
Npars: number of dimensions  
yimp: outputs that need imputed values/outputs with imputed values    
indsp: indices of outputs where output is positive  
indsz: indices of outputs where output is zero  
yp: positive outputs  
Np: # of positive outputs  
xdp: coordinates of positive outputs  
yn: negative/zero outputs  
Nz: # of zero outputs  
xdn: coordinates of negative/zero outputs  
yRL: outputs with imputed values  
indsadd: indices of zero outputs to add to positive points when fitting  
xzs: coordinates of zeros (sorted with the indsadd first)  
xall: all coordinates (sorted with positive points first, then xzs)  
yall: outputs sorted the same way as xall  

RLW_init_impute

xd: design coordinates   
N: # of total points  
indsz: indices of original zero/negative outputs  
xz: coordinates with original z/n outputs  
Nz: # of z/n outputs  
yz: original z/n outputs  
indsp: indices of original positive outputs  
xp: coordinates with original positive outputs  
Np: # of positive outputs  
yp: original z/n outputs  
yc: sample with most negative points  
Npp: # of points that are not zero  
Ncc: # of points that are still zero  
indspp: indices of points that are not zero  
indscc: indices of points that are still zero  
xpp: coordinates of points that aren't zero  
ypp: outputs of xpp  
xcc: coordinates of points that are still zero  
ycc: outputs of xcc  

zGP_gibbs_nrz_optmean

xd: coordinates  
yRL: outputs with imputed values  
Ngibbs: number of iterations to run  
Nzfit: additional # of negative points to include in fit  
Nall: # of total points  
indsz: indices of original zero/negative outputs  
xz: coordinates with original z/n outputs  
Nz: # of z/n outputs  
yz: original z/n outputs  
indsp: indices of original positive outputs  
xp: coordinates with original positive outputs  
Np: # of positive outputs  
yp: original z/n outputs  
ypp: outputs that are positive or negative and included in our fit  
xpp: coordinates of ypp  
xall: coordinates sorted with z/n first  
yall: outputs sorted same as xall   
ysmp: a sample from 1 iteration of outputs  
ysamp: a single output's random sample from inverse CDF method  
ysmpgr: rearranged version of ysmp  
ygr: collection of all samples of the process  
