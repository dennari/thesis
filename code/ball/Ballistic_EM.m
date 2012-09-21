function [opt,lhs,vals] = Ballistic_EM(p0,gi,y,tol_lh,tol_delta,max_iter,min_iter)
global A H m0 u

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
    tol_lh = 1e-4;
end


p = p0;
lh_ = 0;
N = size(y,2)-1; % it is assumed that x0 has y0
lhs = zeros(1,max_iter);
vals = zeros(numel(gi),max_iter);
vals(:,1) = p0(gi)';

% Sigma filter stuff
f = @(x) A*x+repmat(u,1,size(x,2));
h = @(x) H*x;
[usig,w] = CKFPoints(size(A,1));
w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth

for k=1:max_iter
  % E-Step
  %[lh,~,MM,PP,MM_,PP_] = Ballistic_LH(p,y);
  [lh,~,MM,SS,SQ] = Ballistic_LH_Sigma(p,y);
  fprintf(1,'EM %.0f: %.2f %.5f\n',k,lh,exp(p(gi)));
  lhs(k) = lh;
  
  if(k == max_iter); break; end;
  if(min_iter > 0 && k > min_iter)
    if( abs(lh_-lh)             < tol_lh && ... 
        mean(abs(p_(gi)-p(gi))) < tol_delta ); break; end; 
  end
  lh_ = lh;
  p_ = p;
  
  %[MS,PS,DD] = rts_smooth2(MM,PP,MM_,PP_,A);
  %SQ = chol(ballisticQ2D(p(1),p(2)),'lower');
  [MS,SM,DD] = SigmaSmoothSR(MM,SS,f,SQ,usig,w); % D = Smoother Gain
  %[I1,I2,I3] = EM_I123(A,H,m0,y,MS,PS,DD);
  [I1,I2,I3] = EM_I123_Sigma(f,h,m0,y,MS,SM,DD);
  % M-Step
  p = EM_M_Ballistic(p,MS(:,1),gi,N,I1,I2,I3);
  vals(k+1) = p(gi)';
  
  
  
end
lhs = lhs(:,1:k);
vals = vals(:,1:k);
opt = p(gi);








