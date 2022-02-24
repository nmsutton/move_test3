%Ii=0*ones(Ne+Ni,1); % inhibitory input
load('Ii_initial.mat'); % initial gc firing
load('../data/W_Bu09_torus_n900_l2.mat'); % load weight matrix
load('init_firings.mat'); % initial gc firing
mex_hat = W;
t=11;

% generate inhibitory currents
tau = 10; % tau time constant
w_t = zeros(size(mex_hat,1)); % weights multipled by time deltas intermediate values
for i=1:size(Ii)
	% compute weights
	stimes = tbin(i,t,firings);
	for j=1:size(stimes)
        w_t(:,i) = w_t(:,i)+del_t(t-stimes(j));
    end    
end
%w_t = w_t.*mex_hat;
%w_t = (mex_hat*w_t')';
%w_t = w_t.*-0.09;

% calculate tau factor
o = ones(size(mex_hat(:,1)));
w_summed = w_t'*o;
for i=1:size(Ii)
	w_summed2 = 60 - 60*(w_summed/3258);
end
%w_summed = w_summed.*(w_summed<0); % only negative weights
%Ii = w_summed2;
%Ii = Ii.*(Ii>0); % no negative values
w_summed2 = w_summed2.*(w_summed2>0); % no negative values
Ii = Ii + (w_summed2 - Ii)/tau;
Ii_resh = reshape(Ii,30,30);
Ii_resh2 = reshape(w_summed,30,30);

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

function te = del_t(t)
	% value which is reduced according to the time passed (t - spike_time)
	%te = 1/t;
	te = 1/(t^.25);
end