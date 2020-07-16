function [arcoef,ser,resid] = uar(y,n_lags)
resid = NaN(size(y));
% AR coefficients and SER .. no constant term
x = lag(y,1);
for i = 2:n_lags
    x = [x lag(y,i)];
end;
tr_tmp = (1:1:size(y,1))';
tmp = packr([y x tr_tmp]);
itmp = tmp(:,end);
tmp = tmp(:,1:end-1);
y = tmp(:,1);
x = tmp(:,2:end);
bols = inv(x'*x)*(x'*y);
e = y - x*bols;
ssr = e'*e;
ser = sqrt(ssr/(size(x,1)-size(x,2)));
arcoef = bols;
resid(itmp,:) = e;

end

