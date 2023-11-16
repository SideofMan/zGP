function output=zGP_gibbs_nrz_optmean_potential_updates_optimized_clean(xd,yRL,Ngibbs,locs,Nzfit)

locx=locs(1);
locy=locs(2);
y=yRL;

N=length(y);
Nall=N;
Ndim=size(xd,2);

inds=find(y<=0);
xz=xd(inds,:);
indsz=inds;
Nz=length(indsz);
yz=y(indsz,:);

indsp=find(y>0);
Np=length(indsp);
xp=xd(indsp,:);
yp=y(indsp,:);

ypp=y(1:Np+Nzfit);
xpp=xd(1:Np+Nzfit,:);
xall=[xz; xp];
yall=[yz; yp];

sigstarsave=zeros(length(indsz),Ngibbs);

ysmp=yall;
ykeep=ysmp;
yorg=ykeep; % same as y, but order negatives to positives -- goes with xall
ysmpsave=zeros(Nall,Ngibbs);
rhosqsave=zeros(Ngibbs, size(xall,2));

mat52 = @(d) (1+sqrt(5)*d+(5/3)*d.^2).*exp(-sqrt(5)*d); 

ysmp=yall;
if sum(locs)==0
    meantrend='linear_mean';  %here  you could switch up the mean or add a different  mean
else
    meantrend='physical_mean';
end
switch meantrend
    case 'zero_mean'
            %B=mean(ysmp); 
            %mu = @(x) [ones(size(x,1),1)]*B;
            options.zero_mean = true;
            options.nugget_est = false;
            model=ppgasp(xpp,ypp - mean(ypp),options);
    case  'physical_mean'
            % H=@(xd) [ones(size(xd,1),1) xd(:,1) sqrt((xd(:,3)-locx).^2+(xd(:,4)-locy).^2)];
            % 1- vol, 2- BF, 3- E, 4- N -- specific for Aluto
            % 1, vent radius [m]; 2, flux rate (per unit area) [m/s]; 3, bed friction angle [deg]; 4: UTMx of eruptive vent location [m]; 5: UTMy of eruptive vent location [m]
            % Probably  should be
            % H=@(xd) [ones(size(xd,1),1) xd(:,1) xd(:,2) sqrt((xd(:,3)-locx).^2+(xd(:,4)-locy).^2)];
            % B=H(xall)\yall;
            % mu = @(x) H(x)*B;
            options.trend=[ones(size(xpp,1),1) xpp(:,1) xpp(:,2) sqrt((xpp(:,4)-locx).^2+(xpp(:,5)-locy).^2)];
            options.zero_mean = false;
            options.nugget_est = false;
            model=ppgasp(xpp,ypp,options);
    case 'linear_mean'
            options.trend=[ones(Np+Nzfit,1)  xpp];
            options.zero_mean = false;
            options.nugget_est = false;
            model=ppgasp(xpp,ypp,options);
end

% model=ppgasp(xpp,ypp,options);
B=model.theta_hat;
thetas1=model.range_par;
rhosq=thetas1.^2;
rhosq=rhosq';

ykeep=ysmp;
ysmpgr=ysmp;
yorg=ykeep;
count=0;
warning('off','all')

Ball = corr_matrix(xall, xall, rhosq);

Ball=Ball+(1e-6)*eye(size(Ball)); %adding the nugget 1e-3 too big, 1e-12 too small
invBall=inv(Ball);
sigs=(1/Nall)*ysmp'*(Ball\ysmp);
tot_inf_impute=0; %TS
global inf_impute_index; %TSw
inf_impute_index=[(1:Nz)', zeros(Nz,Nz*Ngibbs+1)]; %TS
inf_count = 2; %TS
for ii=1:Ngibbs
    if mod(ii,500)==0
        ii
    end
    % sample range pars every 50th step
    if ii==1 || mod(ii,50)==0
        ii;
        if ii==1
            ypp=y(1:Np+Nzfit);
        else
            ypp=ysmpgr(1:Np+Nzfit);
        end

        model=ppgasp(xpp,ypp,options); %ypps change, so need to find new range pars and B
        B=model.theta_hat; 
        thetas1=model.range_par;
        rhosq=(thetas1.^2)';
        switch meantrend
            case 'zero_mean'
                    B=mean(ysmp); 
                    mu = @(x) [ones(size(x,1),1)]*B;
            case  'physical_mean'
                    H=@(xpp) [ones(size(xpp,1),1) xpp(:,1) xpp(:,2) sqrt((xpp(:,4)-locx).^2+(xpp(:,5)-locy).^2)];
                    % 1- vol, 2- BF, 3- E, 4- N -- specific for Aluto
                    % Probably  should be
                    % H=@(xd) [ones(size(xd,1),1) xd(:,1) xd(:,2) sqrt((xd(:,3)-locx).^2+(xd(:,4)-locy).^2)];
                    %B=H(xall)\yall;
                    mu = @(x) H(x)*B;
            case 'linear_mean'
                    H=@(x)[ones(size(x,1),1),x];
                    mu = @(x) H(x)*B;  
        end

        %%%%%%%%%%%%%%%%%%%%%%%
        Ball = corr_matrix(xpp, xpp, rhosq);
        Ball=Ball+(1e-6)*eye(size(Ball)); %adding the nugget 1e-3 too big, 1e-12 too small
        sigs=(1/(Np+Nzfit))*ypp'*(Ball\ypp);  %calculate sigs using only positive output and close zeros
        %%%%%%%%%%%%%%%%%%%%%% 
        Dsave=zeros(Nall-1,Nz);
        Ballinvsave=zeros(Nall-1,Nall-1,Nz);
        
        Ball = corr_matrix(xall, xall, rhosq);
    
        Ball=Ball+(1e-6)*eye(size(Ball)); %adding the nugget 1e-3 too big, 1e-12 too small
        invBall=inv(Ball);
       
        Ballinv=inv(Ball);
        Ballinvkeep=Ballinv;
        %%%%%%%%%%%%%%%%%%%%%%
        D=ones(Nall-1,1);
        for kk=1:Nz
            if kk==1
                Balltemp=Ball(:,2:end);
                Balltemp=Balltemp(2:end, :);
            elseif kk>1
                Balltemp=Ball(:,[1:kk-1 kk+1:end]);
                Balltemp=Balltemp([1:kk-1 kk+1:end], :);
            end
            if kk==1
                xalltemp=xall(2:end,:);
            elseif kk>1
                xalltemp=xall([1:kk-1 kk+1:end],:);
            end
            for jj=1:Nall-1
                d=sqrt(((xall(kk,:)-xalltemp(jj,:)).^2)./rhosq);
                D(jj)=prod(mat52(d));
            end
            Dsave(:,kk)=D;
            Ballinv=inv(Balltemp);
            Ballinvsave(:,:,kk)=Ballinv;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We'll think about one Gibb's sample as a pass across all 0's where we draw
    % one sample at a time conditioned on all of the other samples. This process
    % ends up will correlated (pass to pass) samples so we will only keep every
    % 5th sample below.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    inf_impute = 0; %TS    
    for jj=1:Nz
        inf_count = inf_count + 1; %TS
        rr=randsample(Nz,1);
        inf_impute_index(rr,5) = inf_impute_index(rr,5)+1;
        ysmp=ykeep;
        if rr==1
            ys=ysmp(2:end);
            xs=xall(2:end,:);
        elseif rr>1
            ys=ysmp([1:rr-1 rr+1:end]);
            xs=xall([1:rr-1 rr+1:end],:);
        end
       
        D=Dsave(:,rr);
        Binv=squeeze(Ballinvsave(:,:,rr));
        sigsqstar=(sigs)*(1-D'*Binv*D);
        sigstarsave(rr,ii)=sigsqstar;
        mustar=mu(xall(rr,:))+D'*(Binv*(ys-mu(xs)));
        % inverse CDF method of sampling a negative value
        yatzero=1/2*(1+erf(-mustar/(sqrt(2)*sigsqstar))); % could add 1e-6 to avoid output being -Inf
        Y=yatzero*rand;
        ypsamp=sqrt(2)*sigsqstar*erfinv(2*Y-1)+mustar;
        if isinf(ypsamp) % need to make sure we have a non infinite value
            %ypsamp = sqrt(2)*sigsqstar*(-5.8636)+mustar;
            ypsamp = -rand/100;
            inf_impute = inf_impute+1; %TS
            inf_impute_index(rr,2)=inf_impute_index(rr,2)+1; %TS
            inf_impute_index(rr,inf_count) = 1; %TS
        end
        yg=ypsamp;
        ykeep(rr)=yg;
    end
    tot_inf_impute=tot_inf_impute+inf_impute;
    sigssave(ii)=sigs;
    ysmp=ykeep;
    ysmpgr(indsz)=ysmp(1:Nz);  % rearranges for fitting range pars
    ysmpgr(indsp)=ysmp(Nz+1:end);
    ysmpsave(:,ii)=ykeep;
    rhosqsave(ii,:)=rhosq;
end
tot_inf_impute %TS
tot_inf_impute/(Nz*Ngibbs) %TS
inf_impute_index(:,1:2) %TS
ygsmps=ysmpsave;
rhogsmps=rhosqsave;

xdr=xd;
xdr(indsz,:)=xall(1:Nz,:);
xdr(indsp,:)=xall(Nz+1:end,:);
ygr=ygsmps;
ygr(indsz,:)=ygsmps(1:Nz,:);
ygr(indsp,:)=ygsmps(Nz+1:end,:);

output{1}=ygr;
output{2}=sqrt(rhogsmps);
output{3}=sigssave;
