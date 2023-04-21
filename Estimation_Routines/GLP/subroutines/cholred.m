function C = cholred(S);
[v,d] = eig((S+S')/2);
d = diag(real(d));
warning off
scale = mean(diag(S))*1e-12;
J = (d>scale);
C = zeros(size(S));
C(J,:) = (v(:,J)*(diag(d(J)))^(1/2))';