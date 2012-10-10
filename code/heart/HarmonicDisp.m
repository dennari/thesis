function [ funs ] = HarmonicDisp()
%HARMONICDISP Summary of this function goes here
%   Detailed explanation goes here

  funs.dispHeader = @dispHeader;
  funs.dispIter = @dispIter;
  funs.dispTrue = @dispTrue;
end

function dispHeader(gi)
  pNames = {'lqw' 'lr' 'lqx'}; pNames = pNames(gi);
  fprintf('Iteration    f(x)        Step-size     %s\n',sprintf('%s      ',pNames{:}));
end

function dispIter(v)
  tmplt = '%5.0f    %13.6g  %13.6g';
  for k=1:numel(v.gi)
    tmplt = [tmplt '  %13.6g'];
  end
  tmplt = [tmplt '\n'];
  
  vv = {v.k,v.lh,v.dx};
  vv(4:4+numel(v.gi)-1) = num2cell(v.p(v.gi));
  
  fprintf(tmplt,vv{:});
  
end




function dispTrue(p_true,gi)
 tmplt = '%5.0f    %13.6g  %13.6g';
  for k=1:numel(gi)
    tmplt = [tmplt '  %13.6g'];
  end
  tmplt = [tmplt '\n'];
  vv = {0,0,0};
  vv(4:4+numel(gi)-1) = num2cell(p_true(gi));
  fprintf('--------\n');
  fprintf(tmplt,vv{:});
end