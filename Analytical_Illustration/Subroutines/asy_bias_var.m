function [bias_var, var_var, var_lp] = asy_bias_var(rho,sigma_2,alpha,h)

    rho2 = rho^2;
    aux = (1+sigma_2^2)/(1-rho2);
    
    bias_var = (h>0).*(alpha.*rho.^(h-1).*(h-1)*sigma_2^2/(aux-1));
    var_var = (h>0).*(rho2.^(h-1).*(1+(h-1).^2/(aux-1))*(1+sigma_2^2)) + rho2.^h*sigma_2^2;
    
    var_lp = aux*(1-rho2.^(h+1))-rho2.^h;

end