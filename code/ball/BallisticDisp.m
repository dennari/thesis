function [ funs ] = BallisticDisp()
%HARMONICDISP Summary of this function goes here
%   Detailed explanation goes here

  funs.dispHeader = @dispHeader;
  funs.dispIter = @dispIter;
  funs.dispTrue = @dispTrue;
  funs.dispInit = @dispInit;
  funs.dispStart = @dispStart;
end

function dispHeader(gi)
  pNames = {'qx' 'qy' 'r' 'qxy' 'ux' 'uy'}; 
  pNames = pNames(gi);
  fprintf('Iteration    f(x)        Step-size     %s\n',sprintf('%s      ',pNames{:}));
end

function dispIter(v)
  tmplt = '%5.0f    %13.6g  %13.6g';
  for k=1:numel(v.gi)
    tmplt = [tmplt '  %13.6g'];
  end
  tmplt = [tmplt '\n'];
  p = v.p; p(3) = exp(p(3));
  vv = {v.k,v.lh,v.dx};
  vv(4:4+numel(v.gi)-1) = num2cell(p(v.gi));
  
  fprintf(tmplt,vv{:});
  
end

function dispInit(p0,gi)
  pNames = {'qx' 'qy' 'r' 'qxy' 'ux' 'uy'};
  p0t = p0; p0t(3) = exp(p0t(3));
  cel = pNames(gi);cel(2,:)=num2cell(p0t(gi));
  fprintf(1,'INITIAL: %s\n',sprintf('%s: %5.4f ',cel{:}));
end

function dispTrue(p_true,gi)
 tmplt = '%5.0f    %13.6g  %13.6g';
  for k=1:numel(gi)
    tmplt = [tmplt '  %13.6g'];
  end
  tmplt = [tmplt '\n'];
  p = p_true; p(3) = exp(p(3));
  vv = {0,0,0};
  vv(4:4+numel(gi)-1) = num2cell(p(gi));
  fprintf('--------\n');
  fprintf(tmplt,vv{:});
end

function dispStart(p_true,p0,gi)
  dispTrue(p_true,gi);
  dispInit(p0,gi);
end




  
