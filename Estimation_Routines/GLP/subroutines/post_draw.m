function [betadraw, sigmadraw] = post_draw(betahat,Sinv,cholZZinv,T)

% Posterior draw of beta and Sigma
% Copied from logMLVAR_formcmc.m in Giannone, Lenza & Primiceri (2015) replication files

n = size(betahat,2);
d=n+2;
eta=mvnrnd(zeros(1,n),Sinv,T+d);
sigmadraw=(eta'*eta)\eye(n);
cholSIGMA=cholred((sigmadraw+sigmadraw')/2);
betadraw=betahat + cholZZinv'*randn(size(betahat))*cholSIGMA;

end

