function [ funs ] = HarmonicDisp()
%HARMONICDISP Summary of this function goes here
%   Detailed explanation goes here

  funs.dispHeader = @dispHeader;
  funs.dispIter = @dispIter;
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
  vv(4:4+numel(v.gi)-1) = num2cell(exp(v.p(v.gi)));
  
  fprintf(tmplt,vv{:});
  
end



