function stimes = tbin(ni, t, firings) 
	% report spike times in a bin of time
	stimes = []; % spike times
	sti = []; % spike time index
	tst = t - 100;% start time of bin
	all_firing = (find(firings(:,2)==ni));
	for si=1:size(all_firing)
		sti = [sti; firings(all_firing(si),1)];
	end
	sti = find(sti>tst & sti<t);
	for si=1:size(sti)
		stimes = [stimes; firings(all_firing(sti(si)),1)];
	end
end