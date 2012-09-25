function [opt,lhs,vals,evals] = Harmonic_BFGS(p0,gi,y,tol_lh,tol_delta,max_iter,min_iter)

if nargin < 7 || isempty(min_iter)
    min_iter = 0;
end
if nargin < 6 || isempty(max_iter)
    max_iter = 1000;
end
if nargin < 5 || isempty(tol_delta)
    tol_delta = 1e-6;
end
if nargin < 4 || isempty(tol_lh)
    tol_lh = 1e-6;
end

opt = optimset(@fminunc);
opt.GradObj = 'on';
opt.LargeScale = 'off'; % use BFGS
opt.TolFun = tol_lh; % force to run for max_iter
opt.TolX = tol_delta;
opt.Display = 'iter'; % prints stuff on every iteration
opt.MaxFunEvals = max_iter;
opt.OutputFcn = @ofun;

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
  [lh,glh] = Harmonic_LH(p,y,gi);
  % remember we're minimizing
  lh = -lh;
  glh = -glh;
end

function stop=ofun(x,val,state)
  %fprintf(1,'BFGS %.0f: %.3f\n',[k,-val.fval]);
  %fprintf(1,'Vals: %.5f\n',exp(x));
  lhs(k) = -val.fval;
  vals(:,k) = x;
  k = k+1;
  stop=0;
end

end









