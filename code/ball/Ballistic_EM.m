function [opt,lhs,vals] = Ballistic_EM(p0,gi,y,tol_lh,tol_delta,max_iter,min_iter)
global A H m0

if nargin < 7 || isempty(min_iter)
    min_iter = 1;
end
if nargin < 6 || isempty(max_iter)
    max_iter = 1000;
end
if nargin < 5 || isempty(tol_delta)
    tol_delta = 1e-3;
end
if nargin < 4 || isempty(tol_lh)
    tol_lh = 1e-6;
end


p = p0;
p_ = p0;
lh_ = 0;
N = size(y,2)-1; % it is assumed that x0 has y0
lhs = zeros(1,max_iter);
vals = zeros(numel(gi),max_iter);
vals(:,1) = p0(gi)';

% Sigma filter stuff
%f = @(x) A*x+repmat(u,1,size(x,2));
%h = @(x) H*x;
%[usig,w] = CKFPoints(size(A,1));
%w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth

for k=1:max_iter
  % E-Step
  [lh,~,MM,PP,MM_,PP_,Q,u] = Ballistic_LH(p,y);
  %[lh,~,MM,SS,SQ] = Ballistic_LH_Sigma(p,y);
%   fprintf(1,'EM %.0f: %.3f\n',[k,lh]);
%   fprintf(1,'Vals: %.5f\n',exp(p(gi)));
  
  lhs(k) = lh;
  lh_delta = lh-lh_;
  x_delta = sqrt(mean((p_(gi)-p(gi)).^2));
  % Display iteration quantities
  fprintf(' %5.0f       %5.0f    %13.6g  %13.6g\n',k,2*k,lh,x_delta);
  
  if(2*k >= max_iter); break; end;
  if(k > 1 && min_iter > 0 && 2*k > min_iter) % don't stop if min_iter not fulfilled
    if( abs(lh_delta)             < tol_lh || ... 
        x_delta                   < tol_delta ); break; end; 
  end
  lh_ = lh;
  p_ = p;
  
  [MS,PS,DD] = rts_smooth2(MM,PP,MM_,PP_,A);
  %SQ = chol(ballisticQ2D(p(1),p(2)),'lower');
  %[MS,SM,DD] = SigmaSmoothSR(MM,SS,f,SQ,usig,w); % D = Smoother Gain
  [I1,I2,I3] = EM_I123(A,H,m0,y,MS,PS,DD,u);
  %[I1,I2,I3] = EM_I123_Sigma(f,h,m0,y,MS,SM,DD);
  % M-Step
  smk = sum(MS(:,2:end),2); smkk = sum(MS(:,1:end-1),2);
  p = EM_M_Ballistic(p,MS(:,1),gi,N,I1,I2,I3,smk,smkk);
  vals(:,k+1) = p(gi)';
  
  
  
end
lhs = lhs(:,1:k);
vals = vals(:,1:k);
opt = p(gi);








