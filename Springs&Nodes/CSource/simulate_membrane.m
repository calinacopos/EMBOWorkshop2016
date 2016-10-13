% 
% Animation of position for motile Dictyostelium cell
% Stepping motility (Mar 2016)
% 

clear;
%close all;
clear all; %clf;

T = load('SideCellN162.txt');
mem = load('membrane.txt');

Npts = 162;
Nframes = length(mem)/Npts;
    
runMode = 'B'; % mode B is animation of position 

figure(1);
for i = 1:Nframes
    % cytoskeleton info
    startPos = (i-1)*Npts + 1; 
    endPos = i*Npts;

    t = mem(startPos:endPos,1);
    x = mem(startPos:endPos,2);
    y = mem(startPos:endPos,3);
    refx = mem(startPos:endPos,4);
    refy = mem(startPos:endPos,5);
    efx = mem(startPos:endPos,6); % elastic force density
    efy = mem(startPos:endPos,7);
    afx = mem(startPos:endPos,8); % adhesive force density
    afy = mem(startPos:endPos,9);
    cfx = mem(startPos:endPos,10);
    cfy = mem(startPos:endPos,11);
    pfx = mem(startPos:endPos,10); % pressre force density
    pfy = mem(startPos:endPos,11);
    fx = mem(startPos:endPos,12); % total force 
    fy = mem(startPos:endPos,13);
    vx = mem(startPos:endPos,18);
    vy = mem(startPos:endPos,19);
    
    switch(runMode)
    case 'B'
        % position figure
        scatter(x,y,'bo','fillcolor','b'); hold on;
        plot(linspace(-15,50,100),-0.5*ones(100,1),'-k');
        xlim([-10 10]);
        ylim([-2 15]);
        set(gca,'plotBoxAspectRatio',[20 17 1]);
        quiver(x,y,fx,fy,'-k');
        quiver(x,y,afx,afy,'-r');
        
        getframe(gcf); %gcf
        hold off;
        %pause 
    end
end;

switch(runMode)
    case 'A'
        figure(2);
        scatter(x,y,'bo','fillcolor','b'); hold on;
        plot(linspace(-15,50,100),-0.5*ones(100,1),'-k');
        hold on;
        axis square;
end