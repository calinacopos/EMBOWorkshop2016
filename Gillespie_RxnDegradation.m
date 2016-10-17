clear;
close all;

K=1; T=50; n=10;
time=zeros(size(T)); 
nn=n*ones(size(T));

for t=1:T % number of events
    T=-(1/(K+nn(t)))*log(rand); 
    time(t+1) = time(t)+T;
    % time for an event, time is updated
    r = rand;
    if r<K/(K+nn(t)) n=n+1; % production
    else n=n-1; % degradation
    end
        nn(t+1)=n; % length is updated
end
figure(1);
stairs(time,nn,'r');