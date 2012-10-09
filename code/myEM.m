function [opt,lhs,vals,times] = myEM(lhf,EM_M,p0,gi,y,dispf,tol_lh,tol_delta,max_iter,min_iter)
global A H f h m0

if nargin < 7 || isempty(min_iter)
    min_iter = 1;
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

p = p0;
p_ = p0;
lh_ = 0;
N = size(y,2)-1; % it is assumed that x0 has y0
lhs = zeros(1,max_iter);
times = lhs;
vals = zeros(numel(gi),max_iter);
vals(:,1) = p0(gi)';

islinear = isempty(f);
[usig,w] = CKFPoints(numel(m0));
w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth



start = tic;
for k=1:max_iter
  %%%%%%%%%%%%%%
  % E-Step%%%%%%
  %%%%%%%%%%%%%%
  if islinear
    [lh,~,MM,PP,MM_,PP_,~,u] = lhf(p,y);
  else
    [lh,~,MM,SS,SQ] = lhf(p,y);
  end

  %%%%%%%%%%%%%%%%%%
  % HOUSEKEEPING %%%
  %%%%%%%%%%%%%%%%%%
  lh_delta = lh-lh_;
  x_delta = sqrt(mean((p_(gi)-p(gi)).^2));
  
  % Display iteration quantities
  dispf(struct('p',p,'gi',gi,'k',k,'lh',lh,'dx',x_delta));

   
  lhs(k) = lh;
  
  if(k >= max_iter); break; end;
  if(k > 1 && min_iter > 0 && k > min_iter) % don't stop if min_iter not fulfilled
    if( abs(lh_delta)             < tol_lh || ... 
        x_delta                   < tol_delta ); break; end; 
  end
  lh_ = lh;
  p_ = p;

  %%%%%%%%%%%%%%%
  % SMOOTH %%%%%%
  %%%%%%%%%%%%%%%
  if islinear
    [MS,PS,DD] = rts_smooth2(MM,PP,MM_,PP_,A);
  else  
    [MS,SM,DD] = SigmaSmoothSR(MM,SS,f,SQ,usig,w); % D = Smoother Gain
  end
  
  try
    if islinear
      [I1,I2,I3] = EM_I123(A,H,m0,y,MS,PS,DD,u);
    else
      [I1,I2,I3] = EM_I123_Sigma(f,h,m0,y,MS,SM,DD);
    end
  catch err
    break;
  end
  
  %%%%%%%%%%%%%%%%
  % M-Step%%%%%%%%
  %%%%%%%%%%%%%%%%
  p = EM_M(p,MS,gi,N,I1,I2,I3);
  
  
  vals(:,k+1) = p(gi)';
  times(k) = toc(start);
  
end
%lhs = lhs(:,1:k);
%vals = vals(:,1:k);
opt = p(gi);








