function tmp=pyplot(o,x,y,s)
	if nargin > 2
		pl.x = x;
		pl.y = y;
		if nargin > 3
			pl.s = s;
		end
	elseif nargin == 2
		pl = x;
	else
		error('at least two arguments required')
  end
  
  if regexp(o,'\.mat$')
    save(o,'-struct','pl','-v7');
  else
    tmp = tempname();
    tmp = [tmp '.mat'];

    save(tmp,'-struct','pl','-v7');
    uns = 'unset LD_LIBRARY_PATH && unset DYLD_LIBRARY_PATH && unset DYLD_FRAMEWORK_PATH && unset OSG_LD_LIBRARY_PATH';
    %pth = '/u/vjvaanan/workspace3/dippa/code';
    python = '/usr/local/bin/python2.7';
    system(sprintf('%s && %s pyplot.py %s %s',uns,python,tmp,o));
  end
end