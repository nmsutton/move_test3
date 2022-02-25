%Ii=0*ones(Ne+Ni,1); % inhibitory input
load('Ii_initial.mat'); % initial gc firing
Ii = Ii_initial;
load('../data/W_Bu09_torus_n900_l2.mat'); % load weight matrix
%load('../data/mex_hat2.mat'); % load weight matrix
load('init_firings.mat'); % initial gc firing
%mex_hat = W;
mex_hat = abs(W);
t=10;

% generate inhibitory currents
tau = 10; % tau time constant
gc_firing = zeros(size(mex_hat,1)); 
for i=1:size(Ii)
	%stimes = tbin(i,t,firings);
    spike_found = find_spike(i,t,firings);
	if spike_found == true
        gc_firing(:,i) = gc_firing(:,i)+1;
    end    
end
in_current = (((mex_hat.^4)*4500)*gc_firing')';

% calculate tau factor
o = ones(size(mex_hat(:,1)));
in_summed = in_current'*o;
%in_summed = gc_firing'*o;
%for i=1:size(Ii)
%	in_summed2 = 150 - 150*(in_summed/900);
%end
in_summed2 = in_summed;
in_summed2 = in_summed2.*(in_summed2>0); % no negative values
Ii = Ii + (in_summed2 - Ii)/tau;
Ii_resh = reshape(Ii,30,30);
Ii_resh2 = reshape(in_summed,30,30);

function spike_found = find_spike(ni, t, firings) 
	% report spike times in a bin of time
	spike_found = false;
	all_spike_times = (find(firings(:,2)==ni));
	for si=1:size(all_spike_times)
		spike_time = firings(all_spike_times(si),1);
        if spike_time == t
            spike_found = true;
        end
    end
end