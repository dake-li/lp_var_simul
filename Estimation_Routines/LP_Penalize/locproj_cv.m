function obj = locproj_cv(obj,lambdaRange)
  
  L = length(lambdaRange);  
  rss = zeros(L,1);
  % aic = zeros(L,1);
  
  for l = 1:L
      % fprintf('.')
      % S = obj.X * inv( obj.X'*obj.X + lambdaRange(l) * obj.P ) * obj.X';
      % rss(l) = sum( ( (obj.Y - S * obj.Y) ./ ( 1 - diag(S) ) ).^2 );
      S_XX = obj.X'*obj.X + lambdaRange(l) * obj.P;
      theta = S_XX \ obj.X' * obj.Y;
      predictY = obj.X * theta;
      diagS = sum(obj.X .* (S_XX \ obj.X')', 2); % leverage
      rss(l) = sum( ( (obj.Y - predictY) ./ ( 1 - diagS ) ).^2 ); % leave-one-out cross validation
  end
  % fprintf('.')
  
  obj.rss     = rss;
  obj.lambdaRange  = lambdaRange;

end