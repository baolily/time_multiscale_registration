
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The linearized Gauss-Seidel relaxation routine for MTV/TV registration
% over the image domain [0,1]x[0,1] or [0,N] x [0,N]
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

function [u1,u2,tau11,tau22,tau12] = mtvgsgm3(N,R,T,u1,u2,g1,g2,IMAX,...
    X1,X2,omega,RiD,Beta,AF)
global alpha_N_fine
global dom
global N_fine

if (dom == 1),
    h = 1/(N);   % the current step size
else
    h = N_fine/N;
end
if (nargin == 14),
    if (RiD==3),
      alpha = AF;  
          D1 = Dmtv(u1,RiD);
          D2 = Dmtv(u2,RiD);      
    else
      alpha = AF/h;
          D1 = Dmtv(u1,RiD);
          D2 = Dmtv(u2,RiD);
    end
else
    if (RiD==3),
     alpha = alpha_N_fine;
         D1 = Dmtv(u1,RiD,Beta);
         D2 = Dmtv(u2,RiD,Beta);
         
%           D1 = myDmtv(u1,u2,RiD,Beta);%% coupling the two components of displacement field 
%           D2 = myDmtv(u2,u1,RiD,Beta);        
    else    
     alpha = alpha_N_fine/h;
         D1 = Dmtv(u1,RiD);
         D2 = Dmtv(u2,RiD);
    end
end
  [n,m] = size(u1);
%       X = [X1(:)';X2(:)'];    
%       U = [u1(:)';u2(:)'];
%     Phi = X+U;    
% X1prime = Phi(1,:)';    
% X2prime = Phi(2,:)';
%      Tk = interp2(X1,X2,T,X1prime,X2prime,'*cubic');
% %      Tk = interp2(X1,X2,T,X1prime,X2prime);
%      Tk = reshape(Tk,size(T));     
% Tk(isnan(Tk)) = 0;

X = [X1(:)';X2(:)'];    
U=[u1(:)';u2(:)'];
Phi=X+U;    X1prime=Phi(1,:)';    X2prime=Phi(2,:)';
Tk = interp2(X1,X2,T,X1prime,X2prime,'*cubic');
Tk = reshape(Tk,size(T));     
Tk(:,1:size(Tk,2)-1) = fillmissing(Tk(:,1:size(Tk,2)-1),'nearest',1);
Tk = fillmissing(Tk,'nearest',2); 

% Tk=movepixels(T,u1,u2);
[TkX1,TkX2] = grad(Tk);
         Dk = R-Tk;  
         f1 = Dk.*TkX1+g1;  
         f2 = Dk.*TkX2+g2;
% if (dom == 1),
      tau11 = TkX1.^2;    
      tau22 = TkX2.^2; 
      tau12 = zeros(N);
% else
%       tau11 = max(TkX1.^2,1.0e-3);    
%       tau22 = max(TkX2.^2,1.0e-3);  
% end         
      oldu1 = u1;     
      oldu2 = u2;
    lambda2 = 1;
    
for iter=1:IMAX %Start the linearized Gauss-Seidel Iteration
    
    u1(1,1)=(1-omega)*u1(1,1)...
        +omega*(f1(1,1)+tau11(1,1)*oldu1(1,1)+alpha*(D1(1,1)*(u1(2,1)+...
        lambda2*u1(1,2))))...
        /(tau11(1,1)+alpha*((1+lambda2)*D1(1,1)));
    
    u2(1,1)=(1-omega)*u2(1,1)...
        +omega*(f2(1,1)+tau22(1,1)*oldu2(1,1)+alpha*(D2(1,1)*(u2(2,1)+...
        lambda2*u2(1,2))))...
        /(tau22(1,1)+alpha*((1+lambda2)*D2(1,1)));
    
    if n>2
    for i=2:n-1
        
    u1(i,1)=(1-omega)*u1(i,1)...
        +omega*(f1(i,1)+tau11(i,1)*oldu1(i,1)+alpha*(D1(i,1)*...
        (u1(i+1,1)+lambda2*u1(i,2))+D1(i-1,1)*u1(i-1,1)))...
        /(tau11(i,1)+alpha*((1+lambda2)*D1(i,1)+D1(i-1,1)));
    
    u2(i,1)=(1-omega)*u2(i,1)...
        +omega*(f2(i,1)+tau22(i,1)*oldu2(i,1)+alpha*(D2(i,1)*...
        (u2(i+1,1)+lambda2*u2(i,2))+D2(i-1,1)*u2(i-1,1)))...
        /(tau22(i,1)+alpha*((1+lambda2)*D2(i,1)+D2(i-1,1)));

    end
    end
    
    u1(n,1)=(1-omega)*u1(n,1)...
        +omega*(f1(n,1)+tau11(n,1)*oldu1(n,1)+alpha*(D1(n,1)*...
        (lambda2*u1(n,2))+D1(n-1,1)*u1(n-1,1)))...
        /(tau11(n,1)+alpha*(lambda2*D1(n,1)+D1(n-1,1)));
    
    u2(n,1)=(1-omega)*u2(n,1)...
        +omega*(f2(n,1)+tau22(n,1)*oldu2(n,1)+alpha*(D2(n,1)*...
        (lambda2*u2(n,2))+D2(n-1,1)*u2(n-1,1)))...
        /(tau22(n,1)+alpha*(lambda2*D2(n,1)+D2(n-1,1)));    
    
    if m>2
    for j=2:m-1   
        
    u1(1,j)=(1-omega)*u1(1,j)...
        +omega*(f1(1,j)+tau11(1,j)*oldu1(1,j)+alpha*(D1(1,j)*...
        (u1(2,j)+lambda2*u1(1,j+1))+D1(1,j-1)*lambda2*u1(1,j-1)))...
        /(tau11(1,j)+alpha*((1+lambda2)*D1(1,j)+lambda2*D1(1,j-1)));
    
    u2(1,j)=(1-omega)*u2(1,j)...
        +omega*(f2(1,j)+tau22(1,j)*oldu2(1,j)+alpha*(D2(1,j)*...
        (u2(2,j)+lambda2*u2(1,j+1))+D2(1,j-1)*lambda2*u2(1,j-1)))...
        /(tau22(1,j)+alpha*((1+lambda2)*D2(1,j)+lambda2*D2(1,j-1)));
    
    for i=2:n-1
        
    u1(i,j)=(1-omega)*u1(i,j)...
        +omega*(f1(i,j)+tau11(i,j)*oldu1(i,j)+alpha*(D1(i,j)*...
        (u1(i+1,j)+lambda2*u1(i,j+1))+D1(i-1,j)*u1(i-1,j)+...
        D1(i,j-1)*lambda2*u1(i,j-1)))...
        /(tau11(i,j)+alpha*((1+lambda2)*D1(i,j)+D1(i-1,j)+...
        lambda2*D1(i,j-1)));
    
    u2(i,j)=(1-omega)*u2(i,j)...
        +omega*(f2(i,j)+tau22(i,j)*oldu2(i,j)+alpha*(D2(i,j)*...
        (u2(i+1,j)+lambda2*u2(i,j+1))+D2(i-1,j)*u2(i-1,j)+...
        D2(i,j-1)*lambda2*u2(i,j-1)))...
        /(tau22(i,j)+alpha*((1+lambda2)*D2(i,j)+D2(i-1,j)+...
        lambda2*D2(i,j-1)));
    
    end 
      
    u1(n,j)=(1-omega)*u1(n,j)...
        +omega*(f1(n,j)+tau11(n,j)*oldu1(n,j)+alpha*(D1(n,j)*...
        (lambda2*u1(n,j+1))+D1(n-1,j)*u1(n-1,j)+D1(n,j-1)*...
        lambda2*u1(n,j-1)))...
        /(tau11(n,j)+alpha*(lambda2*D1(n,j)+D1(n-1,j)+lambda2*D1(n,j-1)));
    
    u2(n,j)=(1-omega)*u2(n,j)...
        +omega*(f2(n,j)+tau22(n,j)*oldu2(n,j)+alpha*(D2(n,j)*...
        (lambda2*u2(n,j+1))+D2(n-1,j)*u2(n-1,j)+D2(n,j-1)*...
        lambda2*u2(n,j-1)))...
        /(tau22(n,j)+alpha*(lambda2*D2(n,j)+D2(n-1,j)+lambda2*D2(n,j-1)));
    end
    end

    u1(1,m)=(1-omega)*u1(1,m)...
        +omega*(f1(1,m)+tau11(1,m)*oldu1(1,m)+alpha*(D1(1,m)*...
        (u1(2,m))+D1(1,m-1)*lambda2*u1(1,m-1)))...
        /(tau11(1,m)+alpha*(D1(1,m)+lambda2*D1(1,m-1)));
    
    u2(1,m)=(1-omega)*u2(1,m)...
        +omega*(f2(1,m)+tau22(1,m)*oldu2(1,m)+alpha*(D2(1,m)*...
        (u2(2,m))+D2(1,m-1)*lambda2*u2(1,m-1)))...
        /(tau22(1,m)+alpha*(D2(1,m)+lambda2*D2(1,m-1)));
    
    if n>2

    for i=2:n-1
        
    u1(i,m)=(1-omega)*u1(i,m)...
        +omega*(f1(i,m)+tau11(i,m)*oldu1(i,m)+alpha*(D1(i,m)*...
        (u1(i+1,m))+D1(i-1,m)*u1(i-1,m)+D1(i,m-1)*lambda2*u1(i,m-1)))...
        /(tau11(i,m)+alpha*(D1(i,m)+D1(i-1,m)+lambda2*D1(i,m-1)));
    
    u2(i,m)=(1-omega)*u2(i,m)...
        +omega*(f2(i,m)+tau22(i,m)*oldu2(i,m)+alpha*(D2(i,m)*...
        (u2(i+1,m))+D2(i-1,m)*u2(i-1,m)+D2(i,m-1)*lambda2*u2(i,m-1)))...
        /(tau22(i,m)+alpha*(D2(i,m)+D2(i-1,m)+lambda2*D2(i,m-1)));
    
    end
    end
    
    u1(n,m)=(1-omega)*u1(n,m)...
        +omega*(f1(n,m)+tau11(n,m)*oldu1(n,m)+alpha*(D1(n-1,m)*...
        u1(n-1,m)+D1(n,m-1)*lambda2*u1(n,m-1)))...
        /(tau11(n,m)+alpha*(D1(n-1,m)+lambda2*D1(n,m-1)));
    
    u2(n,m)=(1-omega)*u2(n,m)...
        +omega*(f2(n,m)+tau22(n,m)*oldu2(n,m)+alpha*(D2(n-1,m)*...
        u2(n-1,m)+D2(n,m-1)*lambda2*u2(n,m-1)))...
        /(tau22(n,m)+alpha*(D2(n-1,m)+lambda2*D2(n,m-1)));
end%End the Gauss-Seidel Iteration