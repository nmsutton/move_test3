% testing a grid cell CAN network with IZ neurons
% Nate Sutton 2022
%Ne=800; Ni=100;
Ne=0; Ni=900;
a=[0.1*ones(Ne,1);0.1+0.08*ones(Ni,1)];
b=[0.2*ones(Ne,1);0.25-0.05*ones(Ni,1)];
c=[-65*ones(Ne+Ni,1)];
d=[8*ones(Ne,1);2*ones(Ni,1)];
p = [a, b, c, d];

v=-65*ones(Ne+Ni,1); % Initial values of v
%Ii=0*ones(Ne+Ni,1); % inhibitory input
load('Ii_initial.mat'); % initial gc firing
Ii = Ii_initial;
u=b.*v;
firings=[];

simdur = 200;%100e3; % total simulation time, ms
ncells = Ne+Ni;%30*30; % total number of cells in network
tau = 10; %% Cell parameters % grid cell synapse time constant, ms
t = 0; % simulation time variable, ms
skip_t = 10; % initial time to skip because pregenerated initial firing is loaded in this time
load('../data/W_Bu09_torus_n900_l2.mat'); % load weight matrix
%load('../data/mex_hat2.mat'); % load weight matrix
load('../data/B_saved.mat'); % velocity input matrix
%load('../data/gc_firing_init.mat'); % initial gc firing
load('init_firings.mat'); % initial gc firing
mex_hat = abs(W);
%Ie=5*B'; % excitatory input
%Ie=60*(B.^2)'; % excitatory input
Ie=60*ones(ncells,1); % excitatory input
%Ii=-1*example_ii';
%Ie = example_ie'*150;
% video parameters
ccol = load('neuron_space_colormap.mat');
savevideo = true;
h = figure('color','w','name','');
numberOfFrames = simdur-skip_t;
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

for t=skip_t:simdur % simulation of 1000 ms
	fired=find(v>=30); % indices of spikes
	firings=[firings; t+0*fired,fired];
	% inhib
	Ii = inhib_curr(Ii, t, mex_hat, firings);
	% excit
	%Ie = Ie .* (1 + (rand(ncells,1)*.02)); % add some random noise
	%Ie = Ie .* (1 - (rand(ncells,1)*.02)); % subtract some random noise
	[v, u] = iznrn(v, u, p, fired, Ie, Ii);
	if savevideo
		myMovie = heatmap(ncells, firings, t, skip_t, h, myMovie, ccol);
	end
end

close(h);
if savevideo
	videofile = VideoWriter('heatmap.avi'); % Create a VideoWriter object to write the video out to a new, different file.
	open(videofile)
	writeVideo(videofile,myMovie) % Write the movie object to a new video file.
	close(videofile)
end

function [v, u] = iznrn(v, u, p, fired, Ie, Ii)
	a=p(:,1);b=p(:,2);c=p(:,3);d=p(:,4);
	v(fired)=c(fired);
	u(fired)=u(fired)+d(fired);
	%v=v+(0.04*v.^2+5*v+140-u+Ie-Ii); % step 1.0 ms
	%v=v+(0.04*v.^2+5*v+140-u+Ie+Ii); % step 1.0 ms
	v=v+(0.04*v.^2+5*v+140-u+Ie-Ii); % step 1.0 ms
	%v=v+(0.04*v.^2+5*v+140-u-2); % step 1.0 ms
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

function Ii = inhib_curr(Ii, t, mex_hat, firings)
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
	in_current = ((mex_hat.^4)*gc_firing')';
    in_current = in_current.*5000;

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
end

function myMovie = heatmap(ncells, firings, t, skip_t, h, myMovie, ccol)
	binned_firing = [];
	for i=1:sqrt(ncells)
		temp = [];
		for j=1:sqrt(ncells)
			nrn_i = ((i-1)*sqrt(ncells))+j;
			spk_t = fbin(nrn_i, t-10, t, firings);
			temp = [temp; spk_t(1,1)];
		end
		binned_firing = [binned_firing; temp'];
	end
	figure(h);
	cla reset;
  	hAxes = gca;
	imagesc(binned_firing);
	colormap(ccol.CustomColormap2);
	xlabel('neuron position on x axis') 
	ylabel('neuron position on y axis')
	shading interp;
	axis square
	title({sprintf('t = %.1f ms',t),'Population activity'})
	set(gca,'ydir','normal')
	cb = colorbar;
	set(cb, 'ylim', [0 5.5]); % set colorbar range
	drawnow
	thisFrame = getframe(gcf);
  	myMovie(t-(skip_t-1)) = thisFrame;
end