function output = probs_zeros(Nz, Np, xdp, yp, xdn, yRL, indsp, indsz)

mat52 = @(d) (1+sqrt(5)*d+(5/3)*d.^2).*exp(-sqrt(5)*d); 
B=zeros(Nz); %correlation matrix of xp points
Btemp=zeros(Nz);

options.trend=[ones(Np,1), xdp];
options.zero_mean = false;
options.nugget_est = false;
modelp=ppgasp(xdp,yp,options);

options.testing_trend=[ones(Nz,1), xdn];
options.mean_only = false;
pred_model=predict_ppgasp(modelp,xdn,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pmean=pred_model.mean;
plot(pmean,'*')
psd=pred_model.sd;
hold on
for k=1:length(pmean)
 line([k k],[pmean(k)-2*psd(k) pmean(k)+2*psd(k)])
 Pn(k)=1/2*(1+erf(-pmean(k)/(psd(k)*sqrt(2))));
end
line([0 Nz],[0 0],'linewidth',3)
  
for j=1:Nz
    for k=1:Np
        distnp(j,k)=sqrt((xdn(j,1)-xdp(k,1)).^2+(xdn(j,2)-xdp(k,2)).^2);
    end
end
mindist=min(distnp, [], 2); % take the minimum across the columns
[valsmd, indsmd]=sort(mindist);
[valspn, indspn]=sort(Pn);
Ntemp=round(.66*length(indspn)); % take the first ~2/3 indices
indsadd=intersect(indspn(1:Ntemp), indsmd(1:Ntemp));

Ninclude=length(indsadd);
indsrest=setdiff(1:1:Nz,indsadd);
xzs=[xdn(indsadd,:); xdn(indsrest,:)];
xall=[xdp; xzs];
yall=[yRL(indsp); yRL(indsz(indsadd)); yRL(indsz(indsrest))];
xp_zp=[xdp; xdn(indsadd,:)];
yp_zp=[yRL(indsp); yRL(indsz(indsadd))];

options.trend=[ones(Np+Ninclude,1)  xp_zp];
model=ppgasp(xp_zp,yp_zp,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
plot(mindist,Pn,'*')
hold on
plot(mindist(indsadd),Pn(indsadd),'r*')

output{1} = xall;
output{2} = yall;
output{3} = Ninclude;