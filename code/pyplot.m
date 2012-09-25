function mat=pyplot(o,data,mat)
  if nargin < 3
    mat = [tempname() '.mat'];
  end
  
  if regexp(o,'\.mat$')
    save(o,'-struct','data','-v7');
  else

    save(mat,'-struct','data','-v7');
    uns = 'unset LD_LIBRARY_PATH && unset DYLD_LIBRARY_PATH && unset DYLD_FRAMEWORK_PATH && unset OSG_LD_LIBRARY_PATH';
    pth = '/Users/dennari/Wrk/dippa/code';
    %python = '/usr/local/bin/python2.7';
    system(sprintf('%s && %s/matpyplot.py %s %s',uns,pth,mat,o));
  end
end