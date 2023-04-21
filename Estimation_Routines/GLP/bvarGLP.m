function r = bvarGLP(y,lags,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function estimates the BVAR of Giannone, Lenza and Primiceri (2012)
%
% y:        data matrix
% lags:     number of lags in the VAR
%
% Last modified: 07/01/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% set BVAR priors (several options available, see setpriors.m)
%  if varargin=[] --> default settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numarg = nargin;
setpriors;


%% data matrix manipulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimensions
[TT,n]=size(y);
k=n*lags+1;         % # coefficients for each equation


% constructing the matrix of regressors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=zeros(TT,k);
% x(:,1)=1;
% for i=1:lags
%     x(:,1+(i-1)*n+1:1+i*n)=lag(y,i);
% end
x=[ones(TT,1) lagmatrix(y,1:lags)];

y0=mean(y(1:lags,:),1);
x=x(lags+1:end,:);
y=y(lags+1:end,:);
[T,n]=size(y);


% MN prior mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=zeros(k,n);
diagb=ones(n,1);
diagb(posi)=0;   % Set to zero the prior mean on the first own lag for variables selected in the vector pos
b(2:n+1,:)=diag(diagb);


%% starting values for the minimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda0=.2;     % std of MN prior
theta0=1;       % std of sur prior
miu0=1;         % std of noc prior
alpha0=2;       % lag-decaying parameter of the MN prior

% residual variance of AR(1) for each variable
SS=zeros(n,1);
for i=1:n
    ar1=ols1(y(2:end,i),[ones(T-1,1),y(1:end-1,i)]);
    SS(i)=ar1.sig2hatols;
end
MIN.psi=SS./100;
MAX.psi=SS.*100;
psi0=SS;

if mn.psi==1;
    inpsi=-log((MAX.psi-psi0)./(psi0-MIN.psi));
elseif mn.psi==0;
    inpsi=[];
end

if mn.alpha==1;
    inalpha=-log((MAX.alpha-alpha0)/(alpha0-MIN.alpha));
elseif mn.alpha==0;
    inalpha=[];
end

if sur==1;
    intheta=-log((MAX.theta-theta0)/(theta0-MIN.theta));
elseif sur==0;
    intheta=[];
end

if noc==1;
    inmiu=-log((MAX.miu-miu0)/(miu0-MIN.miu));
elseif noc==0;
    inmiu=[];
end


x0=[-log((MAX.lambda-lambda0)/(lambda0-MIN.lambda));inpsi;intheta;inmiu;inalpha];

% H0=10*eye(length(x0));          % initial guess for the inverse Hessian



%% maximization of the posterior of the hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [fh,xh,gh,H,itct,fcount,retcodeh] = csminwel('logMLVAR_formin',x0,H0,[],1e-16,1000,y,x,lags,T,n,b,MIN,MAX,SS,Vc,posi,mn,sur,noc,y0,hyperpriors,priorcoef);
fct = @(par) logMLVAR_formin(par,y,x,lags,T,n,b,MIN,MAX,SS,Vc,posi,mn,sur,noc,y0,hyperpriors,priorcoef);
opts = optimoptions('fminunc','Display','notify');
[xh,~,~,output] = fminunc(fct,x0,opts);
itct = output.iterations;


%% output of the maximization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VAR coefficients and residuals at the posterior model
[fh,r.postmax.betahat,r.postmax.sigmahat,r.postmax.T,r.postmax.Sinv,r.postmax.cholZZinv]=fct(xh);

r.lags = lags;                      % # lags
r.postmax.itct=itct;                % #iteration before reaching maximum
r.postmax.SSar1=SS;                 % residual variance of AR(1) for each variable
r.postmax.logPost=-fh;              % value of the posterior of the hyperparameters at the peak
r.postmax.lambda=MIN.lambda+(MAX.lambda-MIN.lambda)/(1+exp(-xh(1)));    % std of MN prior at the peak

r.postmax.theta=MAX.theta;
r.postmax.miu=MAX.miu;

if mn.psi==1;
    % diagonal elements of the scale matrix of the IW prior on the residual variance
    r.postmax.psi=MIN.psi+(MAX.psi-MIN.psi)./(1+exp(-xh(2:n+1)));
    if sur==1;
        % std of sur prior at the peak
        r.postmax.theta=MIN.theta+(MAX.theta-MIN.theta)/(1+exp(-xh(n+2)));
        if noc==1;
            % std of noc prior at the peak
            r.postmax.miu=MIN.miu+(MAX.miu-MIN.miu)/(1+exp(-xh(n+3)));
        end
    elseif sur==0;
        if noc==1;
            % std of sur prior at the peak
            r.postmax.miu=MIN.miu+(MAX.miu-MIN.miu)/(1+exp(-xh(n+2)));
        end
    end
elseif mn.psi==0;
    r.postmax.psi=SS;
    if sur==1;
        % std of sur prior at the peak
        r.postmax.theta=MIN.theta+(MAX.theta-MIN.theta)/(1+exp(-xh(2)));
        if noc==1;
            % std of sur prior at the peak
            r.postmax.miu=MIN.miu+(MAX.miu-MIN.miu)/(1+exp(-xh(3)));
        end
    elseif sur==0;
        if noc==1;
            % std of sur prior at the peak
            r.postmax.miu=MIN.miu+(MAX.miu-MIN.miu)/(1+exp(-xh(2)));
        end
    end
end

if mn.alpha==0;
    r.postmax.alpha=2;
elseif mn.alpha==1;
    % Lag-decaying parameter of the MN prior
    r.postmax.alpha=MIN.alpha+(MAX.alpha-MIN.alpha)/(1+exp(-xh(end)));
end


%% forecasts at the posterior mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Fcast==1
    Y=[y;zeros(hz(end),n)];
    for tau=1:max(hz)
        xT=[1;reshape(Y([T+tau-1:-1:T+tau-lags],:)',k-1,1)]';
        Y(T+tau,:)=xT*r.postmax.betahat;
    end
    r.postmax.forecast=Y(T+hz,:);
end


% %% MCMC
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if mcmc==1
% 
%     % Jacobian of the transformation of the hyperparameters that has been
%     % used for the constrained maximization
%     JJ=exp(xh)./(1+exp(xh)).^2;
%     JJ(1)=(MAX.lambda-MIN.lambda)*JJ(1);
% 
%     if mn.psi==1;
%         JJ(2:n+1)=(MAX.psi-MIN.psi).*JJ(2:n+1);
%         if sur==1;
%             JJ(n+2)=(MAX.theta-MIN.theta)*JJ(n+2);
%             if noc==1;
%                 JJ(n+3)=(MAX.miu-MIN.miu)*JJ(n+3);
%             end
%         elseif sur==0;
%             if noc==1;
%                 JJ(n+2)=(MAX.miu-MIN.miu)*JJ(n+2);
%             end
%         end
%     elseif mn.psi==0;
%         if sur==1;
%             JJ(2)=(MAX.theta-MIN.theta)*JJ(2);
%             if noc==1;
%                 JJ(3)=(MAX.miu-MIN.miu)*JJ(3);
%             end
%         elseif sur==0;
%             if noc==1;
%                 JJ(2)=(MAX.miu-MIN.miu)*JJ(2);
%             end
%         end
%     end
% 
%     if mn.alpha==1;
%         JJ(end)=(MAX.alpha-MIN.alpha)*JJ(end);
%     end
% 
%     JJ=diag(JJ);
%     HH=JJ*H*JJ';
% 
%     % regularizing the Hessian (making sure it is positive definite)
%     [V,E]=eig(HH);
%     HH=V*abs(E)*V';
% 
%     % recovering the posterior mode
%     if mn.psi==1;
%         modepsi=r.postmax.psi;
%     elseif mn.psi==0;
%         modepsi=[];
%     end
% 
%     if mn.alpha==1;
%         modealpha=r.postmax.alpha;
%     elseif mn.alpha==0;
%         modealpha=[];
%     end
% 
%     if sur==1;
%         modetheta=r.postmax.theta;
%     elseif sur==0;
%         modetheta=[];
%     end
% 
%     if noc==1;
%         modemiu=r.postmax.miu;
%     elseif noc==0;
%         modemiu=[];
%     end
% 
%     postmode=[r.postmax.lambda;modepsi;modetheta;modemiu;modealpha];
% 
% 
%     % starting value of the Metropolis algorithm
%     P=zeros(M,length(xh));
%     logMLold=-10e15;
%     while logMLold==-10e15
%         P(1,:)=mvnrnd(postmode,HH*const^2,1);
%         [logMLold,betadrawold,sigmadrawold]=logMLVAR_formcmc(P(1,:)',y,x,lags,T,n,b,MIN,MAX,SS,Vc,pos,mn,sur,noc,y0,...
%             max([MCMCfcast,MCMCstorecoeff]),hyperpriors,priorcoef);
%     end
% 
%     % matrix to store the draws of the VAR coefficients if MCMCstorecoeff is on
%     if MCMCstorecoeff==1
%         r.mcmc.beta  = zeros(k,n,M-N);
%         r.mcmc.sigma = zeros(n,n,M-N);
%     end
% 
%     % matrix to store the forecasts if MCMCfcast is on
%     if MCMCfcast==1
%         r.mcmc.Dforecast=zeros(length(hz),n,M-N);
%     end
% 
%     % Metropolis iterations
%     count=0;
%     for i=2:M
%         if i==100*floor(.01*i);
%                         disp(['Now running the ',num2str(i),'th mcmc iteration (out of ',num2str(M),')'])
%         end
%         % draw candidate value
%         P(i,:)=mvnrnd(P(i-1,:),HH*const^2,1);
%         [logMLnew,betadrawnew,sigmadrawnew]=logMLVAR_formcmc(P(i,:)',y,x,lags,T,n,b,MIN,MAX,SS,Vc,pos,mn,sur,noc,y0,...
%             max([MCMCfcast,MCMCstorecoeff]),hyperpriors,priorcoef);
% 
%         if logMLnew>logMLold
%             logMLold=logMLnew;
%             count=count+1;
%         else
%             if rand(1)<exp(logMLnew-logMLold);
%                 logMLold=logMLnew;
%                 count=count+1;
%             else
%                 P(i,:)=P(i-1,:);
%                 % if MCMCfcast is on, take a new draw of the VAR coefficients
%                 % with the old hyperparameters if have rejected the new ones
%                 % (the speed of this step could be probably improved)
%                 if MCMCfcast==1 | MCMCstorecoeff==1
%                     [junk,betadrawnew,sigmadrawnew]=logMLVAR_formcmc(P(i,:)',y,x,lags,T,n,b,MIN,MAX,SS,Vc,pos,mn,sur,noc,y0,...
%                         max([MCMCfcast,MCMCstorecoeff]),hyperpriors,priorcoef);
%                 end
%             end
%         end
% 
%         % stores draws of VAR coefficients if MCMCstorecoeff is on
%         if i>N & MCMCstorecoeff==1     
%             r.mcmc.beta(:,:,i-N) = betadrawnew;  
%             r.mcmc.sigma(:,:,i-N) = sigmadrawnew;   
%         end
% 
%         % produce and store the forecasts if MCMCfcast is on
%         if i>N & MCMCfcast==1
%             Y=[y;zeros(hz(end),n)];
%             for tau=1:max(hz)
%                 xT=[1;reshape(Y([T+tau-1:-1:T+tau-lags],:)',k-1,1)]';
%                 %Y(T+tau,:)=xT*betadrawnew;
%                 Y(T+tau,:)=xT*betadrawnew+mvnrnd(zeros(1,n),sigmadrawnew);
%             end
%             r.mcmc.Dforecast(:,:,i-N)=Y(T+hz,:);
%             %r.IRF(:,:,:,i) = IRF(betadrawnew,sigmadrawnew);
%         end
% 
%     end
% 
%     % store the draws of the hyperparameters
%     r.mcmc.lambda=P(N+1:end,1);   % std MN prior
% 
%     if mn.psi==1;
%         % diagonal elements of the scale matrix of the IW prior on the residual variance
%         r.mcmc.PSI=P(N+1:end,[2:n+1]);
%         if sur==1;
%             % std of sur prior
%             r.mcmc.theta=P(N+1:end,n+2);
%             if noc==1;
%                 % std of noc prior
%                 r.mcmc.miu=P(N+1:end,n+3);
%             end
%         elseif sur==0;
%             if noc==1;
%                 % std of noc prior
%                 r.mcmc.miu=P(N+1:end,n+2);
%             end
%         end
%     elseif mn.psi==0;
%         if sur==1;
%             % std of sur prior
%             r.mcmc.theta=P(N+1:end,2);
%             if noc==1;
%                 % std of noc prior
%                 r.mcmc.miu=P(N+1:end,3);
%             end
%         elseif sur==0;
%             if noc==1;
%                 % std of noc prior
%                 r.mcmc.miu=P(N+1:end,2);
%             end
%         end
%     end
% 
%     if mn.alpha==1;
%         % Lag-decaying parameter of the MN prior
%         r.mcmc.alpha=P(N+1:end,end);
%     end
% 
%     % acceptance rate
%     r.mcmc.ACCrate=mean(r.mcmc.lambda(2:end)~=r.mcmc.lambda(1:end-1));
% 
% end