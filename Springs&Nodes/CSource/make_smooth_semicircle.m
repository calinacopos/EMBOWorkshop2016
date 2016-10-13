% Make smooth circle
% Calina Copos 
%  
%
R  = 10;    % half width of the base
Nc = 100;
Nb = round(2*Nc/pi);

th = linspace(0,pi,Nc)';
Xc = R*[cos(th),sin(th)];

sb = linspace(-R,R,Nb)';
Xb = [sb(2:Nb-1), 0*sb(2:Nb-1)];

X = [Xb; Xc];

plot(X(:,1),X(:,2),'.-')
hold all;

% diffuse curve to smooth
%
np = size(X,1);
N1x =  floor((np-1)/2);
N2x = (np/2)*ones(rem(np+1,2));
kx =(2/np)* [(0:N1x)  N2x (-N1x:-1)]';

Xhat = fft(X,[],1);
a    = 100;

Xhat = diag(exp( -a*kx.^2 )) * Xhat;

X    = real( ifft( Xhat,[],1));


plot(X(:,1),X(:,2),'.-');
hold off;

% move off 0
% 
X(:,2) = X(:,2)+0.5;

% print to file
%
fid = fopen('SideCellN162.txt','w');
format = '%.16f %.16f \n';
fprintf(fid,format,[X(:,1) X(:,2)]');
fclose(fid);

% The other idea was to use level set curves in electrostatics-like setup
%x = linspace(-2,2,100);
%[x,y] = meshgrid(x);
%g1 = exp(-0.1*(x-0).^2 + (y-0.5).^2);
%g2 = -exp(-0.1*(x-0).^2 + (y+0.5).^2);
%contour(x,y,g1+g2);