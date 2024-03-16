
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A MATLAB code for computing the partial derivative, DIx and DIy of 
% an image I using the central finite dofference formulation with
% Neumann BCs over the image domain [0,1]x[0,1] or [0,N]x[0,N]
%
% DIx and DIy mean the partial derivative of the image I in x- and y
% directions in mathemational sense. In computer programming, we consider 
% the image I as an array I(i,j), where the change in the index i 
% (e.g. i+1 or i-1) means the change in y-direction in mathematical sense. 
% Similarily, the change in the index j (j+1 or j-1) means the change 
% in x-direction in mathematical sense.

% LAST MODIFIED: 2008-August-21
%
% Programed by 
%
% Mr Noppadol Chumchob
% Devision of Applied Mathematics
% Department of Mathematical Sciences
% The University of Liverpool
% Peach Street 
% Liverpool, L69 7ZL, UK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DIx,DIy]=grad(I)
global dom
global N_fine

N = max(size(I)); 
if (dom == 1),
     h = 1/(N);   % the current step size
else
     h = N_fine/N;
end
% [DIx,DIy] = gradient(I,h);
DIx = zeros(N);    DIy = zeros(N);
% Compute the derivatives inside the domain
DIx(2:N-1,2:N-1) = (I(2:N-1,3:N)-I(2:N-1,1:N-2))/(2*h); 
DIy(2:N-1,2:N-1) = (I(3:N,2:N-1)-I(1:N-2,2:N-1))/(2*h);
% Compute the derivatives at the left vertical BC line 
% excluding the corners
DIx(2:N-1,1) = (I(2:N-1,2)-I(2:N-1,1))/(2*h);
DIy(2:N-1,1) = (I(3:N,1)-I(1:N-2,1))/(2*h);
% Compute the derivatives at the right vertical BC line 
% excluding the corners
DIx(2:N-1,N) = (I(2:N-1,N)-I(2:N-1,N-1))/(2*h);
DIy(2:N-1,N) = (I(3:N,N)-I(1:N-2,N))/(2*h);
% Compute the derivatives at the top horizontal BC line 
% excluding the corners
DIx(1,2:N-1) = (I(1,3:N)-I(1,1:N-2))/(2*h);
DIy(1,2:N-1) = (I(2,2:N-1)-I(1,2:N-1))/(2*h);
% Compute the derivatives at the bottom horizontal BC line 
% excluding the corners
DIx(N,2:N-1) = (I(N,3:N)-I(N,1:N-2))/(2*h);
DIy(N,2:N-1) = (I(N,2:N-1)-I(N-1,2:N-1))/(2*h);
% Compute the derivative at the top left
DIx(1,1) = (I(1,2)-I(1,1))/(2*h);
DIy(1,1) = (I(2,1)-I(1,1))/(2*h);
% Compute the derivative at the top right
DIx(1,N) = (I(1,N)-I(1,N-1))/(2*h);
DIy(1,N) = (I(2,N)-I(1,N))/(2*h);
% Compute the derivative at the bottom left
DIx(N,1) = (I(N,2)-I(N,1))/(2*h);
DIy(N,1) = (I(N,1)-I(N-1,1))/(2*h);
% Compute the derivative at the bottom right
DIx(N,N) = (I(N,N)-I(N,N-1))/(2*h);
DIy(N,N) = (I(N,N)-I(N-1,N))/(2*h);
