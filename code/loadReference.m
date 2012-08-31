function [refdata,T] = loadReference(filename)
% LOADREFERENCE - Load physiological reference data from txt files
%
% Syntax:
%   [refdata,T] = loadReference(filename)
%
% In:
%   filename - The physiological data txt file or the folder containing it
%
% Out:
%   refdata - Reference data as a matrix
%   T       - Time points
%
% Description:
% 
%   Load phsyiological reference data acquired during the fMRI scan by
%   respiratory belts etc. This file is specific to the fileformat used by
%   the AMI central devices.
%
%   If the function is called with no output arguments, a figure is shown.
% 
% Example:
%  refpath = '/scratch/ini/datasets/2012-AMI-june-S0340/physiological/00005-EPI-120s';  
%  loadReference(refpath);
%  [refdata,T] = loadReference(refpath);
%  indf = find(refdata(:,6)>0,1,'first');
%  indl = find(refdata(:,6)>0,1,'last');
%  sprintf('Length of scan: %.2f min or s.',T(indl)-T(indf))
%
% Copyright:
%   (c) Arno Solin, 2012
% 


%% Chweck inputs

  % Check argumenents
  if nargin < 1 || isempty(filename),
    error('At least one argument must be specified!')  
  end

  % If the filename is a directory take the first txt file
  if exist(filename,'dir')==7,
      files = dir(fullfile(filename,'*.txt'));
      if length(files) < 1, error('No file found!'); end
      filename = fullfile(filename,files(1).name);
  end
  
  
%% Load and return

  % Load references
  D   = importdata(filename);
  ind = 2:size(D.data,1);
  refdata = D.data(ind,2:end);
  T = D.data(ind,1);
  
  % Show information of the signals
  for i=1:size(refdata,2)
      fprintf('Channel %i was %s: %s given in %s. \n', ...
        i,D.colheaders{1+i},D.textdata{2+2*i,1},D.textdata{3+2*i,1})
  end
  
  
%% Visualize data

if nargout == 0

  figure(1); clf
  for i=1:size(refdata,2);
    subplot(size(refdata,2),1,i)
      plot(T,refdata(:,i))
      title([D.colheaders{1+i} ' ' D.textdata{2+2*i,1}], ...
          'Interpreter','none','FontWeight','bold')
      ylabel(D.textdata{3+2*i,1})
      xlabel('Time [s]')
      ylim auto
  end

end
