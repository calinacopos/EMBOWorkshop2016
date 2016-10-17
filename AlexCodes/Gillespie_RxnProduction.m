clear;
close all;

K=1; T=50; n=10;
time=zeros(size(T)); nn=n*ones(size(T));

for t=1:T % number of events
    tt=-(1/K)*log(rand); 
    time(t+1)=time(t)+tt;
    % time for an event, time is updated
    nn(t+1)=nn(t)+1; % molecule numer is updated
end
figure(2)
stairs(time,nn,'r');

dt = 0.001;
MCtime=(0:dt:T); 
MCnn=n*ones(size(MCtime));
for t=1:T/dt % number of steps
    r=0.5*(1+sign(K*dt-rand)); % random increment for molecule number
    MCnn(t+1)=MCnn(t)+r; % molecule number is updated
end
figure(1)
plot(MCtime,MCnn,'b');