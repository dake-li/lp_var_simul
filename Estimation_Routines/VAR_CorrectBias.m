function [Bc,ByCorrect,Sigma,Sxx] = VAR_CorrectBias(Y,nlags)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nv = size(Y,2); % order (r,x,y,q)
nT = size(Y,1);
[Bc,By,Sigma,Sxx,~] = VAR(Y,nlags);

Ay = reshape(By,[nv,nv * nlags]);
CompAy = kron(diag(ones(1,nlags - 1),-1),eye(nv));
CompAy(1:nv,:) = Ay;
E = eig(CompAy);

sm = zeros(nv * nlags);
for i = 1:(nv * nlags)
    sm = sm + E(i,1) * inv(eye(nv * nlags) - E(i,1) * CompAy');
end

Gamma0 = Sxx(2:end,2:end) - Sxx(2:end,1) * Sxx(2:end,1)';
G = zeros(nv * nlags);
G(1:nv,1:nv) = Sigma;

b = G * (inv(eye(nv * nlags) - CompAy') + CompAy' / (eye(nv * nlags) - CompAy' * CompAy') + sm) / Gamma0;
Bias = - b / (nT - nlags);

CompAyCorrect = CompAy - Bias;
AyCorrect = CompAyCorrect(1:nv,:);
ByCorrect = reshape(AyCorrect,[nv,nv,nlags]);

end

