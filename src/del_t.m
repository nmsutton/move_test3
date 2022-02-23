function te = del_t(t)
	% value which is reduced according to the time passed (t - spike_time)
	%te = 1/t;
	te = 1/(t^.25);
end