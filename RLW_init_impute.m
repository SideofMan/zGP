function output=RLW_init_impute_optimized_clean(xd,y)

warning('off');

N=length(y);

inds=find(y==0);
xz=xd(inds,:);
indsz=inds;
Nz=length(indsz);

inds=find(y>0);
xp=xd(inds,:);
yp=y(inds);
indsp=inds;
Np=length(indsp);
ysave=y;

model=ppgasp(xp,yp-mean(yp)); %j fit a GP to positive values, but make them have mean zero? centering it at zero for simplicity
thetas=(model.range_par);

nugsig=1e-6;
rhosq=((thetas).^2)';

mat52 = @(d) (1+sqrt(5)*d+(5/3)*d.^2).*exp(-sqrt(5)*d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is where we start the negative samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yimputesave=zeros(N,50);
for ii=1:50
    y=ysave;
    if ii==1
        B = corr_matrix(xp, xp, rhosq); % correlation matrix of xp points
        B=B+(nugsig)*eye(size(B)); %adding the nugget (~0) to the diagonal
                                            %of the matrix        
        Btemp = corr_matrix(xd, xd, rhosq);
        Btemp(Np+1:end,:) = 0;
        Btemp(:,Np+1:end) = 0;
        Btemp=Btemp+(nugsig)*eye(size(Btemp)); %adding the nugget (~0) to the diagonal
                                                %of the matrix                                     
        saveB=B;
        saveBtemp=Btemp;
        invB=inv(B);
        saveBinv=invB;
        invBtemp=inv(Btemp);
        saveBtempinv=invBtemp;
        else
           B=saveB;
           Btemp=saveBtemp;
           invB=saveBinv;
           invBtemp=saveBtempinv; 
    end

    r = @(z) mat52(sqrt((repmat(z,[size(xp,1) 1])-xp).^2*(1./rhosq)'));
    ypred = @(z) (r(z))'*(invB*yp);
    sigs=(1/Np)*yp'*(invB*yp);
                     
    Cz = corr_matrix(xz, xz, rhosq); % correlation between censored points (C in doc)
    Czo = corr_matrix(xp, xz, rhosq); % correlation between censored points and positive points (D in doc) 
    
    % better to define ypp here
    ypp=zeros(1,Np);
    for k=1:Np
        % same as before: ypp hasn't been previously defined
        ypp(k) = ypred(xp(k,:));
    end
    mustar=Czo*(invB*(ypp'));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial smps to start RLW alg.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SigSqstar=sigs*(Cz-Czo*(invB*(Czo')));
    Sigma=(SigSqstar+SigSqstar')/2; % assuring Sigma is positive definite
    
    
    smpsize=1e4; % defines the sample size from the multivariate Gaussian
                 % (with #variables equal to Nc)
    smps=mvnrnd(mustar,Sigma,smpsize);
    isneg=zeros(smpsize,Nz);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Batch algorithm for initializing negative-4-zeros sample
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k=1:Nz
        LL=0;
        [indsL, vals]=find(smps(:,k)<LL);
        isneg(indsL,k)=1;
    end
    
    [nsamps,inds]=sort(sum(isneg,2),'descend'); %j sums across columns (outputs rows) and sorts by the samples who had the most J- negatives?
     
    yc=smps(inds(1),:); % makes yc = sample with that most negative points
    inds=find(yc<0);
    y(indsz(inds))=yc(inds); % replace zero outputs with the negative ones from above
    
    
    [inds, vals]=find(y==0); % repeat where we set new indices to the ones that didn't get replaced by the negatives
    indspp=setdiff(1:N,inds);
    Npp=length(indspp);
    Ncc=length(inds);

    indscc=inds;
    xpp=xd(indspp,:);
    xcc=xd(indscc,:);
    ypp=y(indspp);
    ycc=y(indscc);
    count=0;
    while Npp<N % while some of the points are still zero
        count=count+1;
        B = corr_matrix(xpp, xpp, rhosq);
        B=B+(1e-6)*eye(size(B)); % adding the nugget (~0) to the diagonal
                                 % of the matrix 
                                                            
        r = @(z) mat52(sqrt((repmat(z,[size(xpp,1) 1])-xpp).^2*(1./rhosq)'));
        ypred = @(z) (r(z))'*(B\ypp);
        sigs=(1/Npp)*ypp'*(B\ypp);
      
        Cz = corr_matrix(xcc, xcc, rhosq);
        Czo = corr_matrix(xpp, xcc, rhosq);
        
        % better to define ypp here
        ypp=zeros(1,Npp);
        for k=1:Npp
            ypp(k)=ypred(xpp(k,:));
        end
        mustar=Czo*(B\(ypp')); % but ypred seems to include rB^-1(y-mu)...
        SigSqstar=sigs*(Cz-Czo*(B\(Czo')));
        Sigma=(SigSqstar+SigSqstar')/2;
        
        smpsize=10^5; 
        smps=mvnrnd(mustar,Sigma,smpsize);
        isneg=zeros(smpsize,Ncc);
        for k=1:Ncc
            LL=0;
            [indsL,vals]=find(smps(:,k)<LL);
            isneg(indsL,k)=1;
        end
        
        [nsamps,inds]=sort(sum(isneg,2),'descend');
        yctemp=smps(inds(1),:);
        inds=find(yctemp<0);
        y(indscc(inds))=yctemp(inds);
        
        if count>50 % make the remaining zeros very small negatives
            y(indscc)=-rand(size(indscc))/100;
        end

        [inds vals]=find(y==0);
        indspp=setdiff(1:N,inds);
        Npp=length(indspp);
        Ncc=length(inds);
        indscc=inds;
        xpp=xd(indspp,:);
        xcc=xd(indscc,:);
        ypp=y(indspp);
        ycc=y(indscc);
         
    end % end while loop to fill in initial negative values
    yimputestart=y;
    yimputesave(:,ii)=yimputestart;
end
output{1}=yimputesave;
output{2}=sigs;
