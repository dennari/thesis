function pyplot(o,x,y,s)
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

	if isfield(pl,'legend')
		ln = size(pl.legend,2);
		for i=1:ln
			pl.(['legend_' num2str(i)]) = pl.legend{i};
		end
		pl.legend = uint8(ln);
	end

	tmp = tempname();
	save(tmp,'-struct','pl','-v7.3');
	uns = 'unset LD_LIBRARY_PATH && unset DYLD_LIBRARY_PATH && unset DYLD_FRAMEWORK_PATH';
	pth = '/usr/local/share/python';
	system(sprintf('%s && %s/pyplot.py %s.mat %s',uns,pth,tmp,o));
end