#install.packages("R.matlab")
#install.packages("RobustGaSP")
#install.packages("plotly")
#install.packages("MASS")
library("MASS")
library("R.matlab")
library("RobustGaSP")
library("plotly")
library("profvis")

rm(list = ls())
s = Sys.time()
# set wd to xdsave location
set.seed(6)
# toy function #

# xdsave = readMat('BO_toy_design.mat')$xdsave # you could pick these random -- typically that happens with Latin Hypercube sampling, which you can basically thinking of as uniform sampling
# # Toy function adapted from Bastos O'Hagan 2012
# BO_toy = function(x,y) (1-exp(-1/(2*y)))*(2300*x^3+199*x^2+2092*x+60)/(100*x^3+500*x^2+4*x+20)-6
# 
# y = matrix(BO_toy(xdsave[,1], xdsave[,2]))
# xd = xdsave
# ytrue = y

################

# Aluto data # remember to change the locs down low as well

Aluto_data = readMat('Aluto_data_PoIs19_22.mat')

xdsave = Aluto_data$xdtest[[1]][[1]]
y = log10(Aluto_data$ytest[[1]][[1]]+1)
# xdsave = Aluto_data$xdtest[[2]][[1]]
# y = log10(Aluto_data$ytest[[2]][[1]]+1)

maxx = matrix(apply(xdsave, MARGIN = 2, max), nrow = 1)
minx = matrix(apply(xdsave, MARGIN = 2, min), nrow = 1)
diff = maxx - minx

N = dim(xdsave)[1]
Ndim = dim(xdsave)[2]

xd = matrix(0, nrow = N, ncol = Ndim)
for (k in 1:Ndim){
  xd[,k] = (xdsave[,k]-minx[1,k])/diff[1,k]
}

ytrue = y

##################

Nd = dim(xdsave)[1]
Npars = dim(xd)[2]

yimp = y
inds = which(y < 0)
yimp[inds] = 0
ystart = yimp
Ngibbs = 2000

##############################################
# This is where we start the negative samples
##############################################

source("./RLW_init_impute_optimized_clean_fast.R")
source("./corr_matrix.R")

output = RLW_init_impute(xd, ystart) # This does "batch sampling" to get an initially set of negative sampleings
yimputesave = output[[1]]
sigsp = output[[2]]
yimp = matrix(apply(yimputesave, MARGIN = 1, FUN = median))

#########################################################################
# This next  bit is all to arrange order of the design to have [positive
# outputs, closest zeros in design  space]; This is what gets fed into
# zGP_gibbs_nrz_optmean
#########################################################################

y = yimp
N = length(y)
indsp = which(y > 0)
indsz = which(y < 0)
yp = y[indsp,,drop=F]
Np = length(indsp)
xdp = xd[indsp,]
yn = y[indsz,,drop=F]
Nz = length(indsz)
xdn = xd[indsz,]
yRL = matrix(apply(yimputesave, MARGIN = 1, FUN = median))
distnp = matrix(0, nrow = Nz, ncol = Np)

#########################################
# Extra snip to calculate probs of zeros
#########################################

source("./probs_zeros.R")

output = probs_zeros(Nz, Np, xdp, yp, xdn, yRL, indsp, indsz)
xall = output[[1]]
yall = output[[2]]
Ninclude = output[[3]]
ppgasp_options = output[[4]]
erf = output[[5]]
erfinv = output[[6]]

##############################################################################
# Imputing negative responses to design points that have zero outputs via zGP
##############################################################################

source("./zGP_gibbs_nrz_optmean_potential_updates_optimized_clean.R")

# locs = matrix(0, nrow = 1, ncol = 2) # specific for another problem, leave as zeros for now
locs = Aluto_data$locstest[[1]][[1]]
output = zGP_gibbs_nrz_optmean(xall,yall,Ngibbs,locs, Ninclude); # This takes the initial set of negative samples and refines them with Gibbs sampling
                                                                  # output{1} is set of Gibbs samplings for all y
                                                                  # output{2} are (square of) range parameter
                                                                  # samples %j taking original negative samples and resampling with only one point left off

temp = output[[1]]
for (kk in 1:dim(xd)[1]){
  inds[kk] = which(apply(xall, 1, function(row) all(row == xd[kk, ])))
}

#################################################################
# This is the main output of the zGP algorithm. Figure(11) is a
# demonstration of how I imagine it will be used in most cases.
#################################################################

yzgp = matrix(rowMeans(temp[inds,seq(1001,ncol(temp), by = 5)]))

ppgasp_options$trend = cbind(matrix(1, nrow = N, ncol = 1), xd)

modelzgp = ppgasp(design = xd, response = yzgp,
                  trend = ppgasp_options$trend,
                  zero.mean = ppgasp_options$zero.mean,
                  nugget.est = ppgasp_options$nugget.est)
print(Sys.time()-s)
#############################################################
# figure 10: the actual BO_toy function with xdsave plotted
#############################################################

x = seq(0, 1, by = 0.01)
y = seq(0, 1, by = 0.01)
grid = expand.grid(x = x, y = y)
z = with(grid, pmax(BO_toy(grid$x, grid$y), 0))
xx <- matrix(grid$x, nrow = length(x), ncol = length(y), byrow = TRUE)
yy <- matrix(grid$y, nrow = length(x), ncol = length(y), byrow = TRUE)
Ngrid = length(x)
NN = Ngrid*Ngrid

xd = xdsave

plot_ly() %>%
  add_surface(x = x, y = y, z = matrix(z, nrow = length(x), ncol = length(y), byrow = TRUE), alpha = 0.5) %>%
  add_markers(x = xd[,1], y = xd[,2], z = ystart[,])%>%
  layout(title = 'BO_toy surface with xdsave') %>%
  layout(scene = list(aspectmode = "manual", aspectratio = list(x=1,y=1,z=0.5)))

#############################################################
# figure 12: the predicted function with xdsave plotted
#############################################################

yyr = matrix(yy, ncol = 1)
xxr = matrix(xx, ncol = 1)
xyr = cbind(xxr, yyr)

ppgasp_options$testing.trend = cbind(matrix(1, nrow = NN, ncol = 1), xyr)

pred_model = predict.ppgasp(modelzgp, xyr,
                            testing_trend = ppgasp_options$testing.trend)
pmean = pred_model$mean

inds = which(pmean < 0)
pmean[inds] = 0
pmean = matrix(pmean, nrow = Ngrid, ncol = Ngrid)

plot_ly(x = x, y = y, z = pmean, type = "surface", alpha = 0.5) %>%
  add_markers(x = xd[,1], y = xd[,2], z = ystart[,])%>%
  layout(title = 'Predicted surface with xdsave') %>%
  layout(scene = list(aspectmode = "manual", aspectratio = list(x=1,y=1,z=0.5)))

#############################################################
# figure 14: heatmap comparison between the surface and estimate
#############################################################

df = data.frame(xxr, yyr)
df$z = pmax(BO_toy(df$xxr, df$yyr), 0)
df$estimate = pmax(pred_model$mean,0)
df$diff = df$z - df$estimate

x_points = xd[,1]
y_points = xd[,2]
point_df = data.frame(x_points, y_points)

custom_colors = colorRampPalette(colors = c("blue", "turquoise", "orange", "yellow"))(100)

ggplot() + 
  labs(x = "x1", y = "x2", fill = "Difference") +
  geom_tile(data = df, aes(x = xxr, y = yyr, fill = diff)) +
  scale_fill_gradientn(colors = custom_colors) +
  geom_point(data = point_df, aes(x = x_points, y = y_points), color = "black") +
  theme_bw()
