%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The subroutine for computing coeffocients, D, for MTV/TV registration
% over the image domain [0,1]x[0,1] or [0,N]x[0,N]
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
function [D]=myDmtv(u,uu,RiD,Beta)
global beta_N_fine
global dom
global N_fine



 [m,n] = size(u);
if (dom == 1),
     h = 1/(n);   % the current step size
else
     h = N_fine/n;
end
  beta = beta_N_fine*h^2;
    ux = diff(u);
    ux = [ux
zeros(1,m)];
    uy = (diff(u'))';
    uy = [uy zeros(n,1)];
    
    uux = diff(uu);
    uux = [uux
zeros(1,m)];
    uuy = (diff(uu'))';
    uuy = [uuy zeros(n,1)];
    
   if (RiD==3),
%      D = 1./(h^2+Beta*(ux.^2+uy.^2+uux.^2+uuy.^2));    
     D = 1./sqrt((h^2+Beta*(ux.^2+uy.^2+uux.^2+uuy.^2))); 
   else
     D = 1./((sqrt(ux.^2+uy.^2+beta)));
   end
   
   
   
   