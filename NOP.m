%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The routine computing the MTV/TV operator in a dense matrix form  for 
% TV registration over the image domain [0,1]x[0,1] or [0,N]x[0,N]
%
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
function [N] = NOP(u,RiD,Beta,AF)
global alpha_N_fine
global beta_N_fine
global dom
global N_fine

 [m,n] = size(u);

if (dom == 1),
     h = 1/(n);   % the current step size
else
     h = N_fine/n;
end
if (nargin == 4),
    if (RiD==3),
        alpha = AF;
    else
        alpha = AF/h;
    end
else
    if (RiD==3),
        alpha = alpha_N_fine;    
    else
        alpha = alpha_N_fine/h;
    end
end
     beta = beta_N_fine*h^2;
       ux = diff(u);
       ux =([ux
zeros(1,m)]);
       uy = diff(u')';
       uy = ([uy zeros(n,1)]);
    if (RiD==3),
        d = Beta*(ux.^2+uy.^2)+h^2;    
    else
        d = sqrt(ux.^2+uy.^2+beta);
    end
        a = ux./d;
        b = uy./d;
        a = [zeros(1,m)
a];
       ax = (a(2:n+1,:)-a(1:n,:));
        b = [zeros(n,1) b];
       by = (b(:,2:m+1)-b(:,1:m));
        N = -alpha*(ax+by);