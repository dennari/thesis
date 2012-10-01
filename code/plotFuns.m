function [ funs ] = plotFuns()
  funs.normalizeBFGS = @normalizeBFGS;
end

function [lhn,estn] = normalizeBFGS(ev,lh,est)
  %time = time(2:end,:); %remove the initial
  %lh = lh(2:end,:); %remove the initial
  estn = est;
  lhn = lh;
  for j=1:size(lh,2)
    ii = 1;
    df = diff(ev(ev(:,j)>0,j));
    for i=1:numel(df)
      num = df(i);
      if num > 0
        lhn(ii:ii+num-1,j) = interp(lh(i,j),lh(i+1,j),num);
        for kk = 1:size(est,1)
          estn(kk,ii:ii+num-1,j) = interp(est(kk,i,j),est(kk,i+1,j),num); 
        end
        ii = ii + num;
      end
    end
    if ii < size(lh,1)
      if isempty(i)
        i = 0;
      end
      lhn(ii:end,j) = lh(i+1,j); % copy the last one
      for kk = 1:size(est,1)
          estn(kk,ii:end,j) = est(kk,i+1,j); 
      end
    end
  end
end

function l=interp(s,e,kn)
  k=(e-s)/kn;
  l = k*(0:kn-1)+s;
end