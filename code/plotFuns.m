function [ funs ] = plotFuns()
  funs.normalizeBFGS = @normalizeBFGS;
end

function [lh,estn] = normalizeBFGS(time,lh,est,avgt)
  %time = time(2:end,:); %remove the initial
  %lh = lh(2:end,:); %remove the initial
  estn = est;
  for j=1:size(time,2)
    l = lh(:,j);
    ii = 1;
    df = diff(time(time(:,j)>0,j));
    for i=1:numel(df)-1 % discard last one
      num = round(df(i)/avgt);
      if num < 1
        error('num < 1');
      end
      for k=1:num
        l(ii) = lh(i,j); % duplicate lower value num times
        for kk = 1:size(est,1)
          estn(kk,ii,j) = est(kk,i,j); 
        end
        ii = ii + 1;
      end
    end
    l(ii:end) = lh(i+1,j); % copy the last one
    for kk = 1:size(est,1)
        estn(kk,ii:end,j) = est(kk,i+1,j); 
    end
    lh(:,j) = l;
  end
end