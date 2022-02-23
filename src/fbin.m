function spikes = fbin(ni, startt, endt, firings) 
	% report firing in time bin
	% ni = neuron index, startt = start time, endt = end time
	st = []; % spike times
	all_firing = (find(firings(:,2)==ni));
	for si=1:size(all_firing)
		st = [st; firings(all_firing(si),1)];
	end
	st = find(st>startt & st<endt);
	spikes = size(st);
end