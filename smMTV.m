%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A MATLAB code (Demo Version) for computing the smoothing rate for 
% the ralaxation method
%
% LAST MODIFIED: 2008-July-14 
%
% Programed by 
%
% Mr Noppadol Chumchob
% Devision of Applied Mathematics
% Department of Mathematical Sciences
% The University of Liverpool
% Peach Street 
% Liverpool, L69 7ZL, UK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [smooth_rate,low_rate]=smMTV(n,c11,c12,c13,c14,s1,...
    c21,c22,c23,c24,s2,tau11,tau22,tau12,omega)

   n2 = n/2;
   n1 = n/2+1;
   n4 = n/4;
  Low = zeros(n+1,n+1);
 High = zeros(n+1,n+1);
for alpha=-n2:n2
    for beta=-n2:n2
        Eig=RateMTV(alpha,beta,c11,c12,c13,c14,s1,c21,c22,c23,c24,s2,...
            n2,tau11,tau22,tau12,omega);
%         if alpha<=n4 & alpha>=-n4 & beta<=n4 & beta>=-n4
        if alpha<n4 & alpha>-n4 & beta<n4 & beta>-n4            
            High(alpha+n1,beta+n1) = 0;
             Low(alpha+n1,beta+n1) = Eig;
        else
            High(alpha+n1,beta+n1) = Eig;
             Low(alpha+n1,beta+n1) = 0;
        end
    end
end
High = reshape(High,(n+1)*(n+1),1);
 Low = reshape(Low,(n+1)*(n+1),1);
smooth_rate = norm(High,inf);
low_rate = norm(Low,inf);

function [Eig]=RateMTV(alpha,beta,c11,c12,c13,c14,s1,c21,c22,c23,c24,s2,...
    n2,tau11,tau22,tau12,omega) 

AL = zeros(2);

AL(1,1) = (s1+tau11)-omega*(c11*exp(-1/n2*i*pi*alpha)+c13*exp(-1/n2*i*pi*beta));

AL(1,2) = tau12;

AL(2,1) = AL(1,2);

AL(2,2) = (s2+tau22)-omega*(c21*exp(-1/n2*i*pi*alpha)+c23*exp(-1/n2*i*pi*beta));

AR = zeros(2);

AR(1,1) = (1-omega)*(s1+tau11)+omega*(c12*exp(1/n2*i*pi*(alpha))+c14*exp(1/n2*i*pi*(beta)));

AR(1,2) = (1-omega)*tau12;

AR(2,1) = AR(1,2);

AR(2,2) = (1-omega)*(s2+tau22)+omega*(c22*exp(1/n2*i*pi*(alpha))+c24*exp(1/n2*i*pi*(beta)));

A = AL\AR;
EigV = eig(A);

Eig = max(abs(EigV(1)),abs(EigV(2)));