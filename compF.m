%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minor subroutines
%
%   + [f1,f2,Tk,TkX1,TkX2,Dk]=compF(u1,u2,R,T,X1,X2) 
%     <computing force terms>
%
%   + [DIx,DIy]=grad(I)                 <computing image gradient>
%
%   + [fout]=obj_fun(u1,u2,varargin)    <computing the SSD functional>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The subroutine for computing the force term over the image domain 
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
function [f1,f2,Tk,TkX1,TkX2,Dk]=compF(u1,u2,R,T,X1,X2)

% X = [X1(:)';X2(:)'];  
% U = [u1(:)';u2(:)'];
% Phi = X+U;    
% X1prime = Phi(1,:)';    
% X2prime = Phi(2,:)';
% Tk = interp2(X1,X2,T,X1prime,X2prime,'*cubic');
% %邻近点插值Nearest,线性插值Linear,双线性插值bilinear三次样条插值Spline和双立方插值Cubic
% Tk = reshape(Tk,size(T));     Tk(isnan(Tk)) = 0;

X = [X1(:)';X2(:)'];    
U=[u1(:)';u2(:)'];
Phi=X+U;    X1prime=Phi(1,:)';    X2prime=Phi(2,:)';
Tk = interp2(X1,X2,T,X1prime,X2prime,'*cubic');
Tk = reshape(Tk,size(T));     
Tk(:,1:size(Tk,2)-1) = fillmissing(Tk(:,1:size(Tk,2)-1),'nearest',1);
Tk = fillmissing(Tk,'nearest',2); 



% Tk=movepixels(T,u1,u2);
[TkX1,TkX2]=grad(Tk);%modified
Dk = R-Tk; 
f1 = -Dk.*TkX1;    
f2 = -Dk.*TkX2;