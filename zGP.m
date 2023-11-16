%zGP code that goes with Spiller,Wolpert, Tierz, Gasher (2023+) 

clear all
close all
addpath('./functions')
tic;
% toy function

% load BO_toy_design.mat % you could pick these random -- typically that happens with Latin Hypercube sampling, which you can basically thinking of as uniform sampling
% %Toy function adapted from Bastos O'Hagan 2012
% BO_toy = @(x,y) (1-exp(-1./(2*y))).*(2300*x.^3+199*x.^2+2092*x+60)./(100*x.^3+500*x.^2+4*x+20)-6;
% y=BO_toy(xdsave(:,1),xdsave(:,2));
% xd=xdsave;

% ------ %

% Aluto data % remember to change locs down low as well

load Aluto_data_PoIs19_22.mat
xdsave=xdtest{1,1};
y=log10(ytest{1,1}+1);
% xdsave=xdtest{1,2};
% y=log10(ytest{1,2}+1);
xd=xdsave;

% normalizing not needed
% maxx = max(xdsave);
% minx = min(xdsave);
% diff = maxx - minx;
% [N, Ndim] = size(xdsave);
% for k=1:Ndim
%     xd(:,k)=(xdsave(:,k)-minx(k))/diff(k);
% end

% ------ %

ytrue=y;

%%%%%%%%% 
% % testing Aluto data with only 50 points
% 
% testinds = sort(randsample(1:N, 50));
% xd = xd(testinds,:);
% y = y(testinds);

%%%%%%%%%

Nd=size(xd,1);
Npars=size(xd,2);

yimp=y;
inds=find(y<0); % Any x value <0.5 in we'll call 0 -- this should be modified dependent on output under consideration.
yimp(inds)=0;
ystart=yimp;
Ngibbs=2000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is where we start the negative samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output=RLW_init_impute_optimized_clean(xd,ystart); %This does "batch sampling" to get an intitially set of negative sampleings
yimputesave=output{1};
sigsp=output{2};
yimp=median(yimputesave,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This next  bit is all to arrange order of the design to have [positive
%outputs, closest zeros in design  space]; This is what gets fed into
%zGP_gibbs_nrz_optmean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=yimp;
N=length(y);
indsp=find(y>0);
indsz=find(y<0);
yp=y(indsp);
Np=length(indsp);
xdp=xd(indsp,:);
yn=y(indsz);
Nz=length(indsz);
xdn=xd(indsz,:);
yRL=median(yimputesave,2);
distnp=zeros(Nz,Np);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extra snip to calculate probs of zeros
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output = probs_zeros(Nz, Np, xdp, yp, xdn, yRL, indsp, indsz);
xall = output{1};
yall = output{2};
Ninclude = output{3};
global inf_impute_index; %TS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Imputing negative responses to design points that have zero outputs via zGP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%locs=[0 0]; % specific for another problem, leave as zeros for now
locs=locstest{1,1};
output=zGP_gibbs_nrz_optmean_potential_updates_optimized_clean(xall,yall,Ngibbs,locs, Ninclude); %This takes the initial set of negative samples and refines them with Gibbs sampling
                          %output{1} is set of Gibbs samplings for all y
                          %output{2} are (square of) range parameter
                          %samples
                          
temp=output{1};
for kk=1:size(xd,1)  %rearrange response back to original design
    matches = find(all(xd(kk,:)==xall(:,:),2));
    inds(kk)=find(all(xd(kk,:)==xall(:,:),2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main output of the zGP algorithm. Figure(11) is a
% demonstration of how I imagine it will be used in most cases.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yzgp=mean(temp(inds,1001:5:end),2);

options.trend=[ones(N,1), xd];
modelzgp = ppgasp(xd,yzgp,options);
toc
%%
% For plotting
[xx,yy] = meshgrid(0:.01:1, 0:.01:1);
Ngrid=length(xx);
NN=Ngrid*Ngrid;
BOsurf=zeros(101);
xd=xdsave;
for k=1:101
    for j=1:101
        BOsurf(k,j)=max(BO_toy(xx(k,j),yy(k,j)),0);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 10 (the true surface of the toy function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10)
surf(xx,yy,BOsurf,'facealpha',0.5)
shading flat
hold on
plot3(xd(:,1),xd(:,2),ystart,'k*')
xlabel('x_1')
ylabel('x_2')
ah=gca;
set(ah,'fontsize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 12 (the surface obtained by the zGP algorithm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yyr=reshape(yy,NN,1);
xxr=reshape(xx,NN,1);
xyr=[xxr, yyr];

options.testing_trend=[ones(NN,1), xyr];
pred_model=predict_ppgasp(modelzgp,xyr,options);
pmean=pred_model.mean;

pmean=reshape(pmean,Ngrid,Ngrid);
zgpmean=BOsurf;
for k=1:101
    for j=1:101
        zgpmean(k,j)=max(pmean(k,j),0);
    end
end

figure(12)
surf(xx,yy,zgpmean,'facealpha',0.5)
hold on
plot3(xd(:,1),xd(:,2),ystart,'k*')
shading flat
xlabel('x_1')
ylabel('x_2')
ah=gca;
set(ah,'fontsize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 13 (the surface obtained by the zGP algorithm non-censored)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(13)
surf(xx,yy,pmean,'facealpha',0.5)
hold on
plot3(xd(:,1),xd(:,2),yzgp,'k*')
shading flat
xlabel('x_1')
ylabel('x_2')
ah=gca;
set(ah,'fontsize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 14 (the heatmap of the difference between true surface and zGP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(14)
pcolor(xx,yy,BOsurf-zgpmean)
colorbar
shading flat
xlabel('x_1')
ylabel('x_2')
hold on
plot(xd(:,1),xd(:,2),'k*')
shading flat
ah=gca;
set(ah,'fontsize',16)
