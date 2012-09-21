function [opt,lhs,vals,evals] = Ballistic_BFGS(p0,gi,y,tol_lh,tol_delta,max_iter,min_iter)

if nargin < 7 || isempty(min_iter)
    min_iter = 0;
end
if nargin < 6 || isempty(max_iter)
    max_iter = 1000;
end
if nargin < 5 || isempty(tol_delta)
    tol_delta = 1e-2;
end
if nargin < 4 || isempty(tol_lh)
    tol_lh = 1e-6;
end

opt = optimset(@fminunc);
opt.GradObj = 'on';
opt.LargeScale = 'off'; % use BFGS
%opt.TolFun = 1e-30;
%opt.TolX = 1e-30;
opt.Display = 'off';
opt.MaxFunEvals = max_iter;
%opt.OutputFcn = @output;

lhs = zeros(1,max_iter);
vals = zeros(numel(gi),max_iter);
%vals(:,1) = p0(gi)';
k = 1;
[opt,~,~,msg] = fminunc(@bfgs_lh,p0(gi),opt);
lhs = lhs(:,1:k-1);
vals = vals(:,1:k-1);
evals = msg.funcCount; 

function [lh,glh]=bfgs_lh(x)
  p = p0;
  p(gi) = x;
  [lh,glh] = Ballistic_LH(p,y,gi);
  fprintf(1,'BFGS %.0f: %.2f %.5f\n',k,lh,exp(x));
  lhs(k) = lh;
  vals(:,k) = x;
  k = k+1;
  % remember we're minimizing
  lh = -lh;
  glh = -glh;
end


end









