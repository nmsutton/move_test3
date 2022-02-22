% testing a grid cell CAN network with IZ neurons
% Nate Sutton 2022
Ne=800; Ni=100;
a=[0.1*ones(Ne,1);0.1+0.08*ones(Ni,1)];
b=[0.2*ones(Ne,1);0.25-0.05*ones(Ni,1)];
c=[-65*ones(Ne+Ni,1)];
d=[8*ones(Ne,1);2*ones(Ni,1)];
p = [a, b, c, d];

v=-65*ones(Ne+Ni,1); % Initial values of v
Ie=5*ones(Ne+Ni,1); % excitatory input
Ii=0*ones(Ne+Ni,1); % inhibitory input
u=b.*v;
firings=[];

simdur = 200;%100e3; % total simulation time, ms
ncells = Ne+Ni;%30*30; % total number of cells in network
tau = 10; %% Cell parameters % grid cell synapse time constant, ms
t = 0; % simulation time variable, ms
load('../data/W_Bu09_torus_n900_l2.mat'); % load weight matrix
load('../data/B_saved.mat'); % velocity input matrix
load('../data/gc_firing_init.mat'); % initial gc firing
mex_hat = W;
% video parameters
h = figure('color','w','name','');
numberOfFrames = simdur;
allTheFrames = cell(numberOfFrames,1);
vidHeight = 337;%342;
vidWidth = 442;%434;
allTheFrames(:) = {zeros(vidHeight, vidWidth, 3, 'uint8')};
allTheColorMaps = cell(numberOfFrames,1);
allTheColorMaps(:) = {zeros(256, 3)};
myMovie = struct('cdata', allTheFrames, 'colormap', allTheColorMaps);
set(gcf, 'nextplot', 'replacechildren'); 
set(gcf, 'renderer', 'zbuffer');
caxis manual;          % allow subsequent plots to use the same color limits

for t=1:simdur % simulation of 1000 ms
	fired=find(v>=30); % indices of spikes
	firings=[firings; t+0*fired,fired];
	Ii = inhib_curr(Ii, t, mex_hat, firings);
	[v, u] = iznrn(v, u, p, fired, Ie, Ii);
	myMovie = heatmap(ncells, firings, t, h, myMovie);
end
close(h);
%myMovie(1) = []; % remove first frame causing issues due to wrong size
videofile = VideoWriter('heatmap.avi'); % Create a VideoWriter object to write the video out to a new, different file.
open(videofile)
writeVideo(videofile,myMovie) % Write the movie object to a new video file.
close(videofile)

function myMovie = heatmap(ncells, firings, t, h, myMovie)
	binned_firing = [];
	for i=1:sqrt(ncells)
		temp = [];
		for j=1:sqrt(ncells)
			nrn_i = ((i-1)*sqrt(ncells))+j;
			spk_t = fbin(nrn_i, 0, 100, firings);
			temp = [temp; spk_t(1,1)];
		end
		binned_firing = [binned_firing; temp'];
	end
	%custom_colormap = load('neuron_space_colormap.mat');
	figure(h); %ax(1) = subplot(131);
	imagesc(binned_firing);
	colormap(hot);
	axis square
	title({sprintf('t = %.1f ms',t),'Population activity'})
	set(gca,'ydir','normal')
	cb = colorbar;
	set(cb, 'ylim', [0 5.5]); % set colorbar range
	drawnow
	%disp(binned_firing);
	thisFrame = getframe(gcf);
  	myMovie(t) = thisFrame;
end

function [v, u] = iznrn(v, u, p, fired, Ie, Ii)
	a=p(:,1);b=p(:,2);c=p(:,3);d=p(:,4);
	v(fired)=c(fired);
	u(fired)=u(fired)+d(fired);
	v=v+0.5*(0.04*v.^2+5*v+140-u+Ie-Ii); % step 0.5 ms
	v=v+0.5*(0.04*v.^2+5*v+140-u+Ie-Ii); % for numerical stability
	u=u+a.*(b.*v-u);
end

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
	te = 1/t;
end

function Ii = inhib_curr(Ii, t, mex_hat, firings)
	% generate inhibitory currents
	tau = 10; % tau time constant
	w_t = zeros(size(mex_hat,1)); % weights multipled by time deltas intermediate values
	for i=1:size(Ii)
		% compute weights
		w_t = zeros(size(mex_hat,1)); % clear values
		stimes = tbin(i,t,firings);
		for j=1:size(stimes)
			w_t(:,i) = w_t(:,i)+(mex_hat(:,i)*del_t(t-stimes(j))');
		end
	end

	% calculate tau factor
	o = ones(size(mex_hat(:,1)));
	w_summed = w_t*o;
	Ii = Ii + (w_summed - Ii)/tau;
end