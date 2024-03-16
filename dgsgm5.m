%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The coupling pointwise Gauss-Seidel (CPGS) relaxation routine for 
% diffusion registration
% over the image domain [0,1]x[0,1] or [0,N]x[0,N]
%
% LAST MODIFIED: 2008-September-18
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
function [u1,u2,tau11,tau22,tau12] = dgsgm5(N,R,T,u1,u2,g1,g2,IMAX,...
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
      tau12 = TkX1.*TkX2;
      tau21 = tau12;
      tau22 = TkX2.^2;  
% else
%       tau11 = max(TkX1.^2,1.0e-3);    
%       tau22 = max(TkX2.^2,1.0e-3);  
% end         
      oldu1 = u1;           oldu2 = u2;    
         D1 = ones(N);         D2 = ones(N);
       lambda2 = 1;
       
for iter=1:IMAX %Start the Gauss-Seidel Iteration

b=[(f1(1,1)+tau11(1,1)*oldu1(1,1)+tau12(1,1)*oldu2(1,1)...
    +alpha*(D1(1,1)*(u1(2,1)+lambda2*u1(1,2))));...
   (f2(1,1)+tau22(1,1)*oldu2(1,1)+tau21(1,1)*oldu1(1,1)...
   +alpha*(D2(1,1)*(u2(2,1)+lambda2*u2(1,2))))];

a11=(tau11(1,1)+alpha*((1+lambda2)*D1(1,1)));
a12=tau12(1,1);
a21=tau21(1,1);
a22=(tau22(1,1)+alpha*((1+lambda2)*D2(1,1)));

A=[a11 a12;a21 a22];
u=A\b;
u1(1,1)=(1-omega)*u1(1,1)+omega*u(1);
u2(1,1)=(1-omega)*u2(1,1)+omega*u(2);

if n>2

for i=2:n-1

b=[(f1(i,1)+tau11(i,1)*oldu1(i,1)+tau12(i,1)*oldu2(i,1)...
    +alpha*(D1(i,1)*(u1(i+1,1)+lambda2*u1(i,2))+D1(i-1,1)*u1(i-1,1)));...
   (f2(i,1)+tau22(i,1)*oldu2(i,1)+tau21(n,1)*oldu1(i,1)...
   +alpha*(D2(i,1)*(u2(i+1,1)+lambda2*u2(i,2))+D2(i-1,1)*u2(i-1,1)))];

a11=(tau11(i,1)+alpha*((1+lambda2)*D1(i,1)+D1(i-1,1)));
a12=tau12(i,1);
a21=tau21(i,1);
a22=(tau22(i,1)+alpha*((1+lambda2)*D2(i,1)+D2(i-1,1)));

A=[a11 a12;a21 a22];
u=A\b;
u1(i,1)=(1-omega)*u1(i,1)+omega*u(1);
u2(i,1)=(1-omega)*u2(i,1)+omega*u(2);

end

end

b=[(f1(n,1)+tau11(n,1)*oldu1(n,1)+tau12(n,1)*oldu2(n,1)...
    +alpha*(D1(n,1)*(lambda2*u1(n,2))+D1(n-1,1)*u1(n-1,1)));...
   (f2(n,1)+tau22(n,1)*oldu2(n,1)+tau21(n,1)*oldu1(n,1)...
   +alpha*(D2(n,1)*(lambda2*u2(n,2))+D2(n-1,1)*u2(n-1,1)))];

a11=(tau11(n,1)+alpha*(lambda2*D1(n,1)+D1(n-1,1)));
a12=tau12(n,1);
a21=tau21(n,1);
a22=(tau22(n,1)+alpha*(lambda2*D2(n,1)+D2(n-1,1)));

A=[a11 a12;a21 a22];
u=A\b;
u1(n,1)=(1-omega)*u1(n,1)+omega*u(1);
u2(n,1)=(1-omega)*u2(n,1)+omega*u(2);

if m>2

for j=2:m-1

b=[(f1(1,j)+tau11(1,j)*oldu1(1,j)+tau12(1,j)*oldu2(1,j)...
    +alpha*(D1(1,j)*(u1(2,j)+lambda2*u1(1,j+1))+D1(1,j-1)...
    *lambda2*u1(1,j-1)));...
   (f2(1,j)+tau22(1,j)*oldu2(1,j)+tau21(1,j)*oldu1(1,j)...
   +alpha*(D2(1,j)*(u2(2,j)+lambda2*u2(1,j+1))+D2(1,j-1)...
   *lambda2*u2(1,j-1)))];

a11=(tau11(1,j)+alpha*((1+lambda2)*D1(1,j)+lambda2*D1(1,j-1)));
a12=tau12(1,j);
a21=tau21(1,j);
a22=(tau22(1,j)+alpha*((1+lambda2)*D2(1,j)+lambda2*D2(1,j-1)));

A=[a11 a12;a21 a22];
u=A\b;
u1(1,j)=(1-omega)*u1(1,j)+omega*u(1);
u2(1,j)=(1-omega)*u2(1,j)+omega*u(2);

for i=2:n-1
    
b=[(f1(i,j)+tau11(i,j)*oldu1(i,j)+tau12(i,j)*oldu2(i,j)...
    +alpha*(D1(i,j)*(u1(i+1,j)+lambda2*u1(i,j+1))+D1(i-1,j)...
    *u1(i-1,j)+D1(i,j-1)*lambda2*u1(i,j-1)));...
   (f2(i,j)+tau22(i,j)*oldu2(i,j)+tau21(i,j)*oldu1(i,j)...
   +alpha*(D2(i,j)*(u2(i+1,j)+lambda2*u2(i,j+1))+D2(i-1,j)...
   *u2(i-1,j)+D2(i,j-1)*lambda2*u2(i,j-1)))];

a11=(tau11(i,j)+alpha*((1+lambda2)*D1(i,j)+D1(i-1,j)+lambda2*D1(i,j-1)));
a12=tau12(i,j);
a21=tau21(i,j);
a22=(tau22(i,j)+alpha*((1+lambda2)*D2(i,j)+D2(i-1,j)+lambda2*D2(i,j-1)));

A=[a11 a12;a21 a22];
u=A\b;
u1(i,j)=(1-omega)*u1(i,j)+omega*u(1);
u2(i,j)=(1-omega)*u2(i,j)+omega*u(2);
end

b=[(f1(n,j)+tau11(n,j)*oldu1(n,j)+tau12(n,j)*oldu2(n,j)...
    +alpha*(D1(n,j)*(lambda2*u1(n,j+1))+D1(n-1,j)*u1(n-1,j)...
    +D1(n,j-1)*lambda2*u1(n,j-1)));...
   (f2(n,j)+tau22(n,j)*oldu2(n,j)+tau21(n,j)*oldu1(n,j)...
   +alpha*(D2(n,j)*(lambda2*u2(n,j+1))+D2(n-1,j)*u2(n-1,j)...
   +D2(n,j-1)*lambda2*u2(n,j-1)))];

a11=(tau11(n,j)+alpha*(lambda2*D1(n,j)+D1(n-1,j)+lambda2*D1(n,j-1)));
a12=tau12(n,j);
a21=tau21(n,j);
a22=(tau22(n,j)+alpha*(lambda2*D2(n,j)+D2(n-1,j)+lambda2*D2(n,j-1)));

A=[a11 a12;a21 a22];
u=A\b;
u1(n,j)=(1-omega)*u1(n,j)+omega*u(1);
u2(n,j)=(1-omega)*u2(n,j)+omega*u(2);

end

end

b=[(f1(1,m)+tau11(1,m)*oldu1(1,m)+tau12(1,m)*oldu2(1,m)+alpha*(D1(1,m)...
    *(u1(2,m))+D1(1,m-1)*lambda2*u1(1,m-1)));...
   (f2(1,m)+tau22(1,m)*oldu2(1,m)+tau21(1,m)*oldu1(1,m)+alpha*(D2(1,m)...
   *(u2(2,m))+D2(1,m-1)*lambda2*u2(1,m-1)))];

a11=(tau11(1,m)+alpha*(D1(1,m)+lambda2*D1(1,m-1)));
a12=tau12(1,m);
a21=tau21(1,m);
a22=(tau22(1,m)+alpha*(D2(1,m)+lambda2*D2(1,m-1)));

A=[a11 a12;a21 a22];
u=A\b;
u1(1,m)=(1-omega)*u1(1,m)+omega*u(1);
u2(1,m)=(1-omega)*u2(1,m)+omega*u(2);

if n>2

for i=2:n-1
    
b=[(f1(i,m)+tau11(i,m)*oldu1(i,m)+tau12(i,m)*oldu2(i,m)+alpha*(D1(i,m)*...
    (u1(i+1,m))+D1(i-1,m)*u1(i-1,m)+D1(i,m-1)*lambda2*u1(i,m-1)));...
   (f2(i,m)+tau22(i,m)*oldu2(i,m)+tau21(i,m)*oldu1(i,m)+alpha*(D2(i,m)*...
   (u2(i+1,m))+D2(i-1,m)*u2(i-1,m)+D2(i,m-1)*lambda2*u2(i,m-1)))];

a11=(tau11(i,m)+alpha*(D1(i,m)+D1(i-1,m)+lambda2*D1(i,m-1)));
a12=tau12(i,m);
a21=tau21(i,m);
a22=(tau22(i,m)+alpha*(D2(i,m)+D2(i-1,m)+lambda2*D2(i,m-1)));

A=[a11 a12;a21 a22];
u=A\b;
u1(i,m)=(1-omega)*u1(i,m)+omega*u(1);
u2(i,m)=(1-omega)*u2(i,m)+omega*u(2);
end

end

b=[(f1(n,m)+tau11(n,m)*oldu1(n,m)+tau12(n,m)*oldu2(n,m)...
    +alpha*(D1(n-1,m)*u1(n-1,m)+D1(n,m-1)*lambda2*u1(n,m-1)));...
   (f2(n,m)+tau22(n,m)*oldu2(n,m)+tau21(n,m)*oldu1(n,m)...
   +alpha*(D2(n-1,m)*u2(n-1,m)+D2(n,m-1)*lambda2*u2(n,m-1)))];

a11=(tau11(n,m)+alpha*(D1(n-1,m)+lambda2*D1(n,m-1)));
a12=tau12(n,m);
a21=tau21(n,m);
a22=(tau22(n,m)+alpha*(D2(n-1,m)+lambda2*D2(n,m-1)));

A=[a11 a12;a21 a22];
u=A\b;
u1(n,m)=(1-omega)*u1(n,m)+omega*u(1);
u2(n,m)=(1-omega)*u2(n,m)+omega*u(2);

end%End the Gauss-Seidel Iteration

%====================Update tau111 etc
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
      tau12 = TkX1.*TkX2;
      tau21 = tau12;
      tau22 = TkX2.^2; 
%=====================================
