Ne=800; Ni=200;
a=[0.1*ones(Ne,1);0.1+0.08*ones(Ni,1)];
b=[0.2*ones(Ne,1);0.25-0.05*ones(Ni,1)];
c=[-65*ones(Ne+Ni,1)];
d=[8*ones(Ne,1);2*ones(Ni,1)];
p = [a, b, c, d];

v=-65*ones(Ne+Ni,1); % Initial values of v
I=5*ones(Ne+Ni,1); % thalamic input
u=b.*v;
firings=[];
%voltage=[];

for t=1:1000 % simulation of 1000 ms
	fired=find(v>=30); % indices of spikes
	firings=[firings; t+0*fired,fired];
	[v, u] = iznrn(v, u, p, fired, I);
	%voltage(end+1)=v;
end
disp(fire_bin(801,400,500,firings));
%plot(firings(:,1),firings(:,2),'.');

function [v, u] = iznrn(v, u, p, fired, I)
	a=p(:,1);b=p(:,2);c=p(:,3);d=p(:,4);
	v(fired)=c(fired);
	u(fired)=u(fired)+d(fired);
	v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
	v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical stability
	u=u+a.*(b.*v-u);
end

function spikes = fire_bin(ni, startt, endt, firings) 
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