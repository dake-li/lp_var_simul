function obj = locproj_cv(obj,lambdaRange)
  
  L = length(lambdaRange);  
  rss = zeros(L,1);
  % aic = zeros(L,1);
  
  % Preliminary calculations
  S_XX_inv1 = inv(obj.X'*obj.X + obj.P); % Arbitrarily calculate (X'X+lambda*P)^{-1} for lambda=1
  [m,n] = size(obj.D);
  D_XXinv = obj.D*S_XX_inv1(1:n,:);
  D_XXinv_Dp = D_XXinv(:,1:n)*obj.D';
  XY = obj.X'*obj.Y;
  
  for l = 1:L
      % fprintf('.')
      % S = obj.X * inv( obj.X'*obj.X + lambdaRange(l) * obj.P ) * obj.X';
      % rss(l) = sum( ( (obj.Y - S * obj.Y) ./ ( 1 - diag(S) ) ).^2 );
      the_lambda_diff = lambdaRange(l)-1;
      % (X'X+lambda*P)^{-1} via Woodbury formula, noting that P = [D',0]*[D;0]
      S_XX_inv = S_XX_inv1 - the_lambda_diff * D_XXinv' * ( (eye(m) + the_lambda_diff * D_XXinv_Dp) \ D_XXinv );
      theta = S_XX_inv * XY;
      predictY = obj.X * theta;
      diagS = sum(obj.X .* (S_XX_inv * obj.X')', 2); % leverage
      rss(l) = sum( ( (obj.Y - predictY) ./ ( 1 - diagS ) ).^2 ); % leave-one-out cross validation
  end
  % fprintf('.')
  
  obj.rss     = rss;
  obj.lambdaRange  = lambdaRange;

end