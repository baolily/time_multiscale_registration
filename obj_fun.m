
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A MATLAB code for computing the SSD functional over the image domain 
% [0,1]x[0,1] or [0,N]x[0,N]
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
function [fout]=obj_fun(u1,u2,varargin)

R=varargin{1,1}; 
T=varargin{1,2};
X1=varargin{1,3};
X2=varargin{1,4};
      X = [X1(:)';X2(:)'];    
      U = [u1(:)';u2(:)'];
    Phi = X+U;    
X1prime = Phi(1,:)';    
X2prime = Phi(2,:)';
     Tk = interp2(X1,X2,T,X1prime,X2prime,'*cubic');
     Tk = reshape(Tk,size(T));     
Tk(isnan(Tk)) = 0;
tmp=(R-Tk);
r1 = tmp(:);%LEX-Ordering
fout  = 0.5*r1'*r1;
return;