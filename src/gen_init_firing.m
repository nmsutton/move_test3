load('../data/gc_firing_init.mat'); % initial gc firing

bin_size = 10;
max_firing = 5;
max_orig_value = 0.4;
firings = [];

for i = 1:size(gc_firing(1,:)')
	orig_value = gc_firing(i);
	new_val = max_firing * (orig_value/max_orig_value);
	for j = 1:bin_size
		if floor(new_val) > 0
			firing_interval = bin_size/new_val;
			if mod(j,round(firing_interval)) == 0
				temp = [j i];
				firings = [firings; temp];
			end
		end
	end
end