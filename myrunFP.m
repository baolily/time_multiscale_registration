%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A MATLAB code (FP Version 1.0) for registering two given images using 
% the FP methods for the diffusion/MTV/TV registration 
% over the image domain [0,1]x[0,1] or [0,N]x[0,N]
%
% Usage:
%       [u1,u2,Tu,res,t1] = runFP(R,T,RiD,alpha0,omega);
%
% Input parameters:
%
% - R is the reference image (square)
% - T is the template image (square)
% - RiD is the iD for regulariser. 
%       RiD = 1 ---> Diffusion (GS)
%       RiD = 2 ---> Diffusion (CPGS)
%       RiD = 3 ---> MTV       (GS)
%       RiD = 4 ---> TV        (GS)
% - alpha0 >0 is the regularization parameter, usually 1/10 for 
%   the diffusion and MTV models, and 0.25*(1/10) for the TV model. If it 
%   is determined by the cooling approch, but unfortunately the residual 
%   and/or SSD functional are significantly increased during the process, 
%   please use the manual option, setting alpha0 < 1.
% - omega is the relaxation parameter
%   
% Output parameters:
%
% - u1,u2 are 2 components of the deformation u=(u1,u2)
% - Tu is the registered image
% - alpha is the optimal regularization parameter determined by
%     the cooling method. It is designed to give the relative SSD 
%     about 8%-15% for real applications and about 3%-6% for 
%     synthetic examples. 
% - res means the history of the residual
% - t1 is the final runtimes
%
% Ex:
% id = 1; N = 1.0*128; [R,T]=EX(id,N); 
% [u1,u2,Tu,res,t1] = runFP(R,T,1,1/10,1.85);
%
% Remarks :
% 
% - id =1,2
%
% LAST MODIFIED: 2008-September-18 (Final Run for RDR paper)
%
% Deverloped and Programed by 
%
% Noppadol Chumchob and Ke Chen
% Centre for Mathematical Imaging Techniques (CMIT)
% Devision of Applied Mathematics
% Department of Mathematical Sciences
% The University of Liverpool
% Peach Street 
% Liverpool, L69 7ZL, UK    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uh1,uh2,Tk,r2,t1,SSD_V] =myrunFP(R,T,RiD,alpha,omega,Beta,u1,u2)
global alpha_N_fine
global dom
global N_fine 

% Initial process with the given images
[N,M] = size(R);
     
%  Set up the solution system for diffusion registration.
N_fine = N;               %  Size of the original image (Global)
dom = 1;               %  Choice of the domain
                           %  1 --> [0,1]x[0,1], 2 --> [0,N]x[0,N]
GS_iter = 100;             %  No. of the Gauss-Seidel inner iterations 
                           %   (usually 2-3 for real applications)
IMAX_OUT =50;          %  The maximum number of the FAS-cycles
alpha_N_fine = alpha;      %  The global reg. para. 
           
ESP1_OUT = 0.000135;       %  The desired level in tesing the relative 
                           %  reduction of the SSD with respect to the  
                           %  original images, usually 
                           %  1. ESP1_OUT < 0.5 for real applications
                           %      It means that the approach will stop when
                           %      the relative error is less than 50%
                           %  2. ESP1_OUT < 0.05 for syntectic applications
ESP2_OUT = 1.0e-20;        %  The desired level in tesing the convergence 
                           %  of NMG using the change in two consecutive 
                           %  steps of u1 and u2, 1.0e-4 
ESP3_OUT = 1.0e-5;         %  The desired level in tesing the convergence 
                           %  of NGM using the residual obtained from
                           %  2 components of the E-L equation,1.0e-4
ESP4_OUT = 1.0e-10;        %  The desired level in tesing the convergence 
                           %  of NGM using the change in two consecutive 
                           %  steps (or Grad of) of the fitting term 
                           %  (default = 1.0e-5 [dom=1], 0.1 [dom=2])
ESP5_OUT = 1.0e-20;        %  The desired level in tesing the convergence 
                           %  of NGM using the change in two consecutive 
                           %  steps of the relative registered image
                           %   (default = 1.0e-2)
done_out = false;          %  The outer loop indicator
k = 1;                     %  The counter for the outer iteration
w = [0.05 0.25 0.40 0.25 0.05]; 
W = w'*w;                  % Gaussian filter for making smooth images
                           % before starting registration process
SM = 1;                    % The option for computing the smoothing rate 
                           % during deformation process, 
                           % SM = 0 -->Compute, SM = 1 --> Not Compute
ACC_SMV = ones(N);                           

% Original input images                           
 R0 = R; T0 = T;     
 
% The smooth version of R and T
%   R = conv2(R,W,'same');    T = conv2(T,W,'same');

if (dom == 1),
    h_fine = 1/(N_fine);    % The step size for discretization the image domain
                        % [0,1]x[0,1] based on the cell-centred approach
    x1d = ((1:N_fine)-1)*h_fine; %Positions in the discreted domain in
    x2d = x1d;              % X1- and X2- directions 
else
    h_fine = 1;             % The step size for discretization the image domain
                        % [0,N]x[0,N] based on the cell-centred approach
    x1d = (2.*(1:N_fine)-1)*(h_fine/2); %Positions in the discreted domain in
    x2d = x1d;              % X1- and X2- directions 
end

[X1,X2]= meshgrid(x1d,x2d); % Grid positions over the discreted domain 
X = [X1(:)';X2(:)'];   

% Initial solutions
% - uh1 and uh2 is the initial guess solution for the first - and second-
%   displacment field u=(u1,u2)
    uh1 = zeros(N);    uh2 = zeros(N);
%     uh1 = u1;    uh2 = u2;
    gh1 = zeros(N);    gh2 = zeros(N);
% Computing the partial derivative of the new transformed template image
    Tk = T;
    [TkX1,TkX2] = grad(Tk);
    Dk = Tk-R ;
% Computing the new force
    fh1 = -Dk.*TkX1;  fh2 = -Dk.*TkX2;  
% Computing the initial residual
    if ((RiD == 1) | (RiD == 2)),
        R1 = gh1-(LOP(uh1)+fh1);  
        R2 = gh2-(LOP(uh2)+fh2);    
    elseif (RiD == 3),
        R1 = gh1-(NOP(uh1,RiD,Beta)+fh1);    %NOP--laplacian operation  
        R2 = gh2-(NOP(uh2,RiD,Beta)+fh2);    
    elseif (RiD == 4),
        R1 = gh1-(NOP(uh1,RiD,Beta)+fh1);    
        R2 = gh2-(NOP(uh2,RiD,Beta)+fh2);    
    end
 res10 = norm(R1(:),2);
 res20 = norm(R2(:),2);  
 r2(k) = max(res10,res20);
% Compute the SSD between two input images
% - SSD0 is the initial SSD obatined from two input images R and T
% - SSD is the SSD values
    %SSD0  = (h_fine^2)*obj_fun(zeros(N),zeros(N),R,T,X1,X2); %The first SSD 1/2*(R-T)^2
    SSD0 = 1/2*(norm(R(:)-T(:),2))^2;
    SSD   = SSD0;
     SSD_V(1)=SSD;
%     SSD_V(k)  = SSD; 
    NT0   = norm((T(:)),2);
    NT(k) = NT0;  %The first relative Tk   
%  Perform the non-linear PDE solver.
 %   fprintf('---------FP for the diffusion/MTV/TV registration--------\n');
 %   fprintf('The number of grids at the finest level is %d x %d \n',N,N);
    fprintf('---------------------------------------------------------\n');
    fprintf('Initial SSD = %g \n',SSD0);    
  %  fprintf('Initial ||N1(u1)-F1(U)|| = %g \n',res10);    
  %  fprintf('Initial ||N2(u2)-F2(U)|| = %g \n',res20);    
  %  fprintf('Initial ||N(U)-F(U)|| = %g \n\n',max(res10,res20));
    tic;
while(~done_out),
    OLDuh1=uh1;   %OLDuh1 is the previous first component of the 
                  %   displacement field u
    OLDuh2=uh2;   %OLDuh2 is the previous second component of the 
                  %   displacement field u
    oldTk=Tk;     %oldTk is the previous deformed template image
% Call the FP routine
    if (RiD == 1),
% w-GS         
    [uh1,uh2,tau11,tau22,tau12] = dgsgm3(N,R,T,uh1,uh2,gh1,gh2,GS_iter,...
        X1,X2,omega);
    elseif (RiD == 2),
% w-CPGS         
    [uh1,uh2,tau11,tau22,tau12] = dgsgm5(N,R,T,uh1,uh2,gh1,gh2,GS_iter,...
        X1,X2,omega);
    elseif (RiD == 3),
% w-GS        
    [uh1,uh2,tau11,tau22,tau12] = mtvgsgm3(N,R,T,uh1,uh2,gh1,gh2,...
        GS_iter,X1,X2,omega,RiD,Beta);
% w-CPGS    
%     [uh1,uh2,tau11,tau22,tau12] = mtvgsgm5(N,R,T,uh1,uh2,gh1,gh2,...
%         GS_iter,X1,X2,omega,RiD,Beta);
    elseif (RiD == 4),
    [uh1,uh2] = mtvgsgm3(N,R,T,uh1,uh2,gh1,gh2,GS_iter,X1,X2,omega,RiD);
    end
% Computing the new transformed template image
% Applying the MATLAB commands for the interpolation process in wraping the
% template image and computing the new force
% [fh1,fh2,Tk,TkX1,TkX2]=compF(uh1,uh2,R,T,X1,X2);
    [fh1,fh2,Tk] = compF(uh1,uh2,R,T,X1,X2);
% Compute the residual for equations
if ((RiD == 1) | (RiD == 2)),
   r1m = gh1-(LOP(uh1)+fh1);  
   r2m = gh2-(LOP(uh2)+fh2); 
elseif (RiD == 3),
   r1m = gh1-(NOP(uh1,RiD,Beta)+fh1);  
   r2m = gh2-(NOP(uh2,RiD,Beta)+fh2); 
elseif (RiD == 4),
   r1m = gh1-(NOP(uh1,RiD)+fh1);  
   r2m = gh2-(NOP(uh2,RiD)+fh2); 
end
if (dom == 1),
    res1 = norm(r1m(:),2)/(res10);
    res2 = norm(r2m(:),2)/(res20);
else
    res1 = norm(r1m(:),2);
    res2 = norm(r2m(:),2);
end      
% Computing the new SSD
   %SSD = (h_fine^2)*obj_fun(uh1,uh2,R,T,X1,X2);
%    SSD_V(1)=SSD;
   SSD=1/2*(norm(R(:)-Tk(:),2))^2;
   
     k = k+1;
% - res1out and res2out are the difference between solutions
   res1out(k) = norm(OLDuh1(:)-uh1(:),2);% res1out--sqrt((oldu1-u1).^2)
   res2out(k) = norm(OLDuh2(:)-uh2(:),2);
% - res3out and res4out are the residual between the left- and right-hand
%    sides
   res3out(k) = res1;
   res4out(k) = res2;
% -r1 and r2 are the maximum of the differences between solutions and 
%   the residuals
% -J is the current SSD
    r1(k) = max(res1out(k),res2out(k));
    r2(k) = max(res3out(k),res4out(k));
    SSD_V(k) = SSD;  
    NT(k) = norm(oldTk(:)-Tk(:),2); 
    
    
% fprintf('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n')
% fprintf('Iteration = %d--SSD = %f, Err = %f  \n',k-1,...
%     SSD,SSD/SSD0)%Err--the relative  of sum of squared difference
% fprintf('||OLDu1 - NEWu1|| = %g \n',res1out(k));
% fprintf('||OLDu2 - NEWu2|| = %g \n',res2out(k));
% %fprintf('||N1(u1)-F1(U)||  = %g \n',res1);
% %fprintf('||N2(u2)-F2(U)||  = %g \n\n',res2);
% %fprintf('||OLDsol(U)-NEWsol(U)|| = %g \n',r1(k));
% %fprintf('||N(U)-F(U)|| = %g \n',r2(k));
% fprintf('||OLDssd - NEWssd|| = %g \n',abs(SSD_V(k)-SSD_V(k-1)));
% fprintf('||OLD_Tk - NEW_Tk|| = %g \n',abs(NT(k)-NT(k-1)));
% fprintf('CPU - Time =  %g s\n\n',toc);
% fprintf('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n')
% % Computing the smoothing rate

if SM == 0,
    if (RiD>3),
        fprintf('The smoother is not the Gauss-Seidel type...\n')
        t1=toc;
        time=t1
        return
    end
  SM_Rate = zeros(N);
%Diffusion Setting  
if (RiD==1),
    c11 = alpha;  c21 = alpha;
    c12 = alpha;  c22 = alpha;
    c13 = alpha;  c23 = alpha;
    c14 = alpha;  c24 = alpha;
     s1 = c11+c12+c13+c14;   
     s2 = c21+c22+c23+c24;
  if (dom==1),
    h_fine = 1/N;
    c11 = alpha/h_fine^2;  c21 = alpha/h_fine^2;
    c12 = alpha/h_fine^2;  c22 = alpha/h_fine^2;
    c13 = alpha/h_fine^2;  c23 = alpha/h_fine^2;
    c14 = alpha/h_fine^2;  c24 = alpha/h_fine^2;
     s1 = c11+c12+c13+c14;   
     s2 = c21+c22+c23+c24;
  end
end
%MTV Setting
if (RiD>2),
%    D1 = Dmtv(uh1,RiD,Beta);
%    D2 = Dmtv(uh2,RiD,Beta); 
 
      D1 = myDmtv(uh1,uh2,RiD,Beta);%% coupling the two components of displacement field 
      D2 = myDmtv(uh2,uh1,RiD,Beta); 
   
end
for i=2:N-1,
    for j=2:N-1
% Smoothing Rate MTV
if (RiD>2),
c11 = alpha*D1(i-1,j);  c21 = alpha*D2(i-1,j);
c12 = alpha*D1(i+1,j);  c22 = alpha*D2(i+1,j);
c13 = alpha*D1(i,j-1);  c23 = alpha*D2(i,j-1);
c14 = alpha*D1(i,j+1);  c24 = alpha*D2(i,j+1);
 s1 = c11+c12+c13+c14;   s2 = c21+c22+c23+c24;
end
        SM_Rate(i,j) = smMTV(N,c11,c12,c13,c14,s1,c21,c22,c23,c24,s2,...
            tau11(i,j),tau22(i,j),tau12(i,j),omega);        
    end
end
            SMR(k-1) = max(max(SM_Rate(2:N-1,2:N-1)));
ACC_SMV(2:N-1,2:N-1) = ACC_SMV(2:N-1,2:N-1)...
    .*((SM_Rate(2:N-1,2:N-1)).^GS_iter);
         ACC_SM(k-1) = max(max(ACC_SMV(2:N-1,2:N-1)));  
        ACC_SMR(k-1) = max(max(SM_Rate(2:N-1,2:N-1).^GS_iter)); 
    
fprintf('SMu  = %g, ACC.SMR = %g, ACC_SMV = %g \n',SMR(k-1),ACC_SMR(k-1),...
    ACC_SM(k-1));
end

% Checking the convergence for solution
if (SSD/SSD0<ESP1_OUT) 
    fprintf('\nThe algorithm has reached the desired relative ');
    fprintf('error at %g percent in %d iterations\n',(SSD/SSD0)*100,k-1);
    t1=toc;
    time=t1
    done_out=true;
elseif (k>1 & r1(k)<ESP2_OUT)
    t1=toc;
    cpu_time = t1
    fprintf('\nThe solution converges!(The SOL. does not change any more) ');
    fprintf('at level %g\n',ESP2_OUT);
    done_out=true;
elseif (k>=2 & r2(k)<ESP3_OUT)
    t1=toc;
    cpu_time = t1
    fprintf('\nThe avg. norm for ||N(U)-F(U)|| has been reached the ');
    fprintf('desired level! at level %g \n',ESP3_OUT);
    done_out=true;
elseif (k>1 & abs(SSD_V(k)-SSD_V(k-1))<ESP4_OUT)
    t1=toc;
    cpu_time = t1
    fprintf('\nThe SSD convergesd at the level = %g\n',ESP4_OUT);
    done_out=true;
elseif (k>1 & abs(NT(k)-NT(k-1))<ESP5_OUT)
    t1=toc;
    cpu_time = t1
    fprintf('\nThe registered image converges at level = %g\n',....
        ESP5_OUT);
    done_out=true;
% elseif (k>1 & SSD_V(k)-SSD_V(k-1)>0)
%     t1=toc;
%     cpu_time = t1
%     fprintf('\nstoped k=%d \n',....
%         k);
%     done_out=true;
    
elseif (k>IMAX_OUT)
    t1=toc;
    cpu_time = t1
    fprintf('\n The maximum number iteration is exceed %d \n',IMAX_OUT);
    done_out=true;
end
% keyboard
end %while loop
% Visualize the all results at the end
%    close all
%     figure;
%     subplot( 1,3,1);
%     colormap gray;
%     imagesc( R0 );
%     title('Reference image R')
%     xlabel(['\alpha = ',num2str(alpha)])
%     axis on;
%     axis square;

%     subplot( 1,3,2);
%     colormap gray;
%     imagesc( T0 );
%     axis on;
%     axis square;
%     title('Template image T')
%     xlabel(['SSD = ',num2str(SSD0),', Error = ',...
%         num2str(100*(SSD0/SSD0)),'%'])

%     U=[uh1(:)';uh2(:)'];
%     Phi=X+U;    X1prime=Phi(1,:)';    X2prime=Phi(2,:)';
%     Tk = interp2(X1,X2,T0,X1prime,X2prime,'*cubic');
%     Tk = reshape(Tk,size(T));     
%     Tk(:,1:size(Tk,2)-1) = fillmissing(Tk(:,1:size(Tk,2)-1),'nearest',1);
%     Tk = fillmissing(Tk,'nearest',2); 
    SSD_V=SSD_V(1:end-1);
    
%     Tk(isnan(Tk)) = 0;   
%     subplot( 1,3,3);
%     colormap gray;
%     imagesc(Tk); 
%     title(['Registered Image T_{k}']);
% %     title(['Registered Image T_{k}, k = ',num2str(k-1),...
% %         ', Time = ',num2str(t1/60),' mins']);
%     xlabel(['SSD = ',num2str(SSD),', Error =',...
%         num2str(100*(SSD/SSD0)),'%'])
%     axis on;
%     axis square;

%% plot SSD VS Iteration
%     figure
% %     subplot(2,2,1)
%     semilogy(1:k,SSD_V,'r-.d','LineWidth',2)
%     xlim([1 k]);
%     title('SSD VS. No. Iteration');
%     xlabel(['Number of Iterations      (Time =',num2str(toc/60),' mins)']);
%     %ylabel('D[U]');
%      ylabel('SSD');


%%

%     subplot(2,2,2)
%     semilogy([0 1:k-1],[1 res3out(2:k)],...
%        'r-s',[0 1:k-1],[1 res4out(2:k)],'b-.d','LineWidth',2)
%     h = legend('1st eq.','2nd eq.',1);
%     if (k>2)
%     xlim([1 k-1]);
%     end
%     set(h,'Interpreter','none')
%     title('||OLDu - NEWu|| and ||N(U) - F(U)|| VS. No. Iteration');
%     xlabel(['Number of Iterations    (Time = ',num2str(toc/60),' mins )']);
%     ylabel('||OLDu - NEWu|| and ||N(U) - F(U)||');
%     
%     subplot(2,2,3)
%     mesh(abs(r1m))
%     title('Absolute value of the residual for 1st eq.');    
% 
%     subplot(2,2,4)
%     mesh(abs(r2m))
%     title('Absolute value of residual for 2nd eq.');    




%     u=sqrt(uh1.^2+uh2.^2);
%     figure;
%     colormap gray;
%     imagesc( u );
%     title('Magnitude of the deformation field u');
%     axis square
%     colorbar
%%%%%%%%%%%%%%%%%%%%%%%%% END MAIN CONTROL PROGRAM%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Main Subroutines
%
% + [u1,u2] = dgsgm3(N,R,T,u1,u2,g1,g2,IMAX,X1,X2,omega,AF)
%   <GS for diffusion registration> 
%
% + [u1,u2] = dgsgm5(N,R,T,u1,u2,g1,g2,IMAX,X1,X2,omega,AF)
%   <CPGS for diffusion registration> 
%
% + [L] = LOP(u,AF)
%   <computing the dense matrix for Laplacian>
%
% + [u1,u2] = mtvgsgm3(N,R,T,u1,u2,g1,g2,IMAX,X1,X2,omega,RiD,AF)
%   <GS for MTV/TV registration>
%
% + [D]=Dmtv(u,RiD)
%   <computing the coefficients for MTV/TV >
% 
% + [N] = NOP(u,RiD,AF)
%   <computing the dense matrix for MTV/TV operator>
% 













