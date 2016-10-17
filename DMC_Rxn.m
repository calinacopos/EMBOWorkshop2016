% Direct Monte-Carlo simulation for reaction equation
% Presented by Alex Mogilner (Oct 17, 2016)
%
%

K=100; dt=0.001; tend=8;
t=(0:dt:tend); P=zeros(size(t));
for s=1:(tend/dt)
    P(s+1) = P(s)+0.5*(1+sign(K*dt-rand))-0.5*(1+sign(P(s)*dt-rand));
end

mean(P)
std(P)

figure(1);
plot(t,P);

figure(2); 
hist(P,30)