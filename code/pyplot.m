function pyplot(o,x,y,s)
	pl.x = x;
	pl.y = y;
	if nargin > 3
		pl.s = s;
	end
	tmp = tempname();
	save(tmp,'-struct','pl','-v7.3');
	uns = 'unset LD_LIBRARY_PATH && unset DYLD_LIBRARY_PATH && unset DYLD_FRAMEWORK_PATH';
	pth = '/usr/local/share/python';
	system(sprintf('%s && %s/pyplot.py %s.mat %s',uns,pth,tmp,o));
end