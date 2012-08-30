function plotstruct(s,ax)
	if nargin < 2
		ax = gca();
	end
	for k=1:size(s.data,2)
		d = s.data{k};
		if size(d,2) > 2
			plot(ax,d{1},d{2},d{3});
		else
			plot(ax,d{1},d{2});
		end
		hold on;
	end
	if isfield(s,'xlabel')
		xlabel(ax,s.xlabel);
	end
	if isfield(s,'ylabel')
		ylabel(ax,s.ylabel);
	end
	if isfield(s,'title')
		title(ax,s.title);
	end
	if isfield(s,'legend')
		legend(ax,s.legend);
	end
end