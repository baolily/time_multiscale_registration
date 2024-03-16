%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The linearized Gauss-Seidel relaxation routine for diffusion registration
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
function [u1,u2,tau11,tau22,tau12] = dgsgm3(N,R,T,u1,u2,g1,g2,IMAX,...
    X1,X2,omega,AF)
global alpha_N_fine
global dom
global N_fine

if (dom == 1),
    h = 1/(N);   % the current step size
else
    h = N_fine/N;
end
if (nargin == 12), 
  alpha = AF/h^2;
else
  alpha = alpha_N_fine/h^2;
 end
  [n,m] = size(u1);
      X = [X1(:)';X2(:)'];    
      U = [u1(:)';u2(:)'];
    Phi = X+U;    
X1prime = Phi(1,:)';    
X2prime = Phi(2,:)';
     Tk = interp2(X1,X2,T,X1prime,X2prime,'*cubic');
     Tk = reshape(Tk,size(T));     
Tk(isnan(Tk)) = 0;
[TkX1,TkX2] = grad(Tk);
         Dk = R-Tk;  
         f1 = Dk.*TkX1+g1;  
         f2 = Dk.*TkX2+g2;
% if (dom == 1),
      tau11 = TkX1.^2;
      tau12 = zeros(N);
      tau22 = TkX2.^2;  
% else
%       tau11 = max(TkX1.^2,1.0e-3);    
%       tau22 = max(TkX2.^2,1.0e-3);  
% end         
      oldu1 = u1;     
      oldu2 = u2;          

for iter=1:IMAX %Start the linearized Gauss-Seidel Iteration
    
    u1(1,1)=(1-omega)*u1(1,1)+omega/(2*alpha+tau11(1,1))*(f1(1,1)...
        +tau11(1,1)*oldu1(1,1)+alpha*(u1(2,1)+u1(1,2)));
    
    u2(1,1)=(1-omega)*u2(1,1)+omega/(2*alpha+tau22(1,1))*(f2(1,1)...
        +tau22(1,1)*oldu2(1,1)+alpha*(u2(2,1)+u2(1,2)));

    if n>2
    for i=2:n-1
        
    u1(i,1)=(1-omega)*u1(i,1)+omega/(3*alpha+tau11(i,1))*(f1(i,1)...
        +tau11(i,1)*oldu1(i,1)+alpha*(u1(i+1,1)+u1(i,2)+u1(i-1,1)));
    
    u2(i,1)=(1-omega)*u2(i,1)+omega/(3*alpha+tau22(i,1))*(f2(i,1)...
        +tau22(i,1)*oldu2(i,1)+alpha*(u2(i+1,1)+u2(i,2)+u2(i-1,1)));        
    
    end
    end
    
    u1(n,1)=(1-omega)*u1(n,1)+omega/(2*alpha+tau11(n,1))*(f1(n,1)...
        +tau11(n,1)*oldu1(n,1)+alpha*(u1(n,2)+u1(n-1,1)));
    
    u2(n,1)=(1-omega)*u2(n,1)+omega/(2*alpha+tau22(n,1))*(f2(n,1)...
        +tau22(n,1)*oldu2(n,1)+alpha*(u2(n,2)+u2(n-1,1)));
    
    if m>2
    for j=2:m-1   
        
    u1(1,j)=(1-omega)*u1(1,j)+omega/(3*alpha+tau11(1,j))*(f1(1,j)...
        +tau11(1,j)*oldu1(1,j)+alpha*(u1(2,j)+u1(1,j+1)+u1(1,j-1)));
    
    u2(1,j)=(1-omega)*u2(1,j)+omega/(3*alpha+tau22(1,j))*(f2(1,j)...
        +tau22(1,j)*oldu2(1,j)+alpha*(u2(2,j)+u2(1,j+1)+u2(1,j-1)));  
    
    for i=2:n-1
        
    u1(i,j)=(1-omega)*u1(i,j)+omega/(4*alpha+tau11(i,j))*(f1(i,j)...
        +tau11(i,j)*oldu1(i,j)+alpha*(u1(i+1,j)+u1(i,j+1)+u1(i-1,j)...
        +u1(i,j-1)));
    
    u2(i,j)=(1-omega)*u2(i,j)+omega/(4*alpha+tau22(i,j))*(f2(i,j)...
        +tau22(i,j)*oldu2(i,j)+alpha*(u2(i+1,j)+u2(i,j+1)+u2(i-1,j)...
        +u2(i,j-1)));
    
    end 
      
    u1(n,j)=(1-omega)*u1(n,j)+omega/(3*alpha+tau11(n,j))*(f1(n,j)...
        +tau11(n,j)*oldu1(n,j)+alpha*(u1(n,j+1)+u1(n-1,j)+u1(n,j-1)));
    
    u2(n,j)=(1-omega)*u2(n,j)+omega/(3*alpha+tau22(n,j))*(f2(n,j)...
        +tau22(n,j)*oldu2(n,j)+alpha*(u2(n,j+1)+u2(n-1,j)+u2(n,j-1)));
    end
    end

    u1(1,m)=(1-omega)*u1(1,m)+omega/(2*alpha+tau11(1,m))*(f1(1,m)...
        +tau11(1,m)*oldu1(1,m)+alpha*(u1(2,m)+u1(1,m-1)));
    
    u2(1,m)=(1-omega)*u2(1,m)+omega/(2*alpha+tau22(1,m))*(f2(1,m)...
        +tau22(1,m)*oldu2(1,m)+alpha*(u2(2,m)+u2(1,m-1)));
    
    if n>2

    for i=2:n-1
        
    u1(i,m)=(1-omega)*u1(i,m)+omega/(3*alpha+tau11(i,m))*(f1(i,m)...
        +tau11(i,m)*oldu1(i,m)+alpha*(u1(i+1,m)+u1(i-1,m)+u1(i,m-1)));
    
    u2(i,m)=(1-omega)*u2(i,m)+omega/(3*alpha+tau22(i,m))*(f2(i,m)...
        +tau22(i,m)*oldu2(i,m)+alpha*(u2(i+1,m)+u2(i-1,m)+u2(i,m-1)));       
    
    end
    end
    
    u1(n,m)=(1-omega)*u1(n,m)+omega/(2*alpha+tau11(n,m))*(f1(n,m)...
        +tau11(n,m)*oldu1(n,m)+alpha*(u1(n-1,m)+u1(n,m-1)));
    
    u2(n,m)=(1-omega)*u2(n,m)+omega/(2*alpha+tau22(n,m))*(f2(n,m)...
        +tau22(n,m)*oldu2(n,m)+alpha*(u2(n-1,m)+u2(n,m-1)));
end%End the Gauss-Seidel Iteration