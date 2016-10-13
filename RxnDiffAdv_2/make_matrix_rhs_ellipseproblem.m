%
% MAKE_MATRIX_RHS_ELLIPSEPROBLEM -- make the matrix corresponding to the 
%                                  discrete Laplacian on the unit square 
%                                  excluding a specified circle inside
%                                  the square f is the rhs corresponding
%                                  to u=1 on the circle equations are
%                                  solved everywhere, inside and out,
%                                  for simplicity
%                   
%     input:  n   -- number of interior grid points in each direction
%
%     output: A   -- (n*n)x(n*n) matrix of modified discrete Laplacian
%             f   -- (n*n)x1 right-hand-side vector
%             phi -- signed distance function, phi<0 inside circle
%                                              phi>0 outside circle
%
function [A,f,phi]=make_matrix_rhs_ellipseproblem(n)
    
    % initialize A as the standard Laplacian (scaled)
    %
    e = ones(n,1);
    L1 = spdiags([e -2*e e],-1:1,n,n);
    I1 = speye(n);
    A  = kron(L1,I1) + kron(I1,L1);
    
    % initalize f to be zeros
    %
    f = zeros(n*n,1);

    % boundary value assumed to be constant
    %
    %Ub=0;
    %Ub=0.5;
    
    % form arrays of grid point locations
    %
    dx = 2/(n-1);
    x = linspace(-1,1,n); 
    [x,y] = ndgrid(x);
   
    % parameters for the embedded ellipse
    %
    xc = 0.0;
    yc = 0.0;
    rad = 1.0;
    a = 0.2;
    b = 0.6;

    % compute the signed distance function
    %
    phi = sqrt( ((x-xc).^2)./a^2 + ((y-yc).^2)./b^2 ) - rad;
   
    % loop over the interior points and flag the irregular points
    %
    IJ = [-1  0;
           1  0;
           0 -1;
           0  1;
         ];
    
    for j=2:n-1
        for i=2:n-1
            
            % skip the interior points
            %
            if( phi(i,j) > 0 )
                continue
            end
                
            % check the neighbors
            %ai
            for k=1:4
                
               if( phi(i+IJ(k,1),j+IJ(k,2)) > 0 )
                   
                   % approximate distance to the boundary, scaled by h
                   %
                   alpha = phi(i,j)/( phi(i,j) - phi(i+IJ(k,1),j+IJ(k,2)));
                   
                   % compute the distance to the boundary
                   %
                   kr = sub2ind([n,n],i,j);
                   kc = sub2ind([n,n],i+IJ(k,1),j+IJ(k,2));
                   
                   %keyboard
                   
                   % adjust RHS
                   %
                   %Ub = y(i,j); % BC is u = sin(theta) = y (for r=1)
                   Ub = 0.25*(1+25*(1-0.32*y(i,j)*y(i,j)));
                   f(kr) = f(kr) - Ub/alpha;
                   
                   % adjust diagonal element
                   %
                   A(kr,kr) = A(kr,kr) + 1 - 1/alpha;
                   
                   % adjust off-diagional, and enforce symmetry
                   %
                   A(kr,kc) = 0;
                   A(kc,kr) = 0;
                   
               end
            end
               
        end
    end
    