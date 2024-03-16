%% The energy function is updated by template image and alpha & beta  every time
%   update 2020-02-04 
clc
close all
clear 

%% Reading the two images %%%%%
R = imread('synthetic_reference_32.png');
T = imread('synthetic_template_32.png');

%% Initializing 
N=size(T,1);
R = double(R); 
MaxR = 1/max(max(R)); R = R*MaxR;
T = double(T); 
MaxT = 1/max(max(T)); T = T*MaxT; 

x1d = (1:N)./N; 
x2d = x1d; 
[X1,X2]= meshgrid(x1d,x2d);
 X = [X1(:)';X2(:)']; 
 
Omega=1.85;
Tk=T;
ux=zeros(N);uy=zeros(N);
SSD=[];

%% Solving
for kk=1:7
     beta=5000*kk;
     alpha=100*0.5^kk; 

    [u1,u2,Tk,res,t1,SSD_V] = myrunFP(R,Tk,3,alpha,Omega,beta);
    Alpha(kk)=alpha;
    Beta(kk)=beta;
    ux=ux+u1;
    uy=uy+u2;
    SSD=[SSD SSD_V];
end

%% Printing the result before and after registration
SSD_before =1/2*(norm(R(:)-T(:),2))^2;
SSD_after = 1/2*(norm(R(:)-Tk(:),2))^2;
fprintf('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n')
fprintf('SSD_before = %f--SSD_after = %f, Err = %f  \n',SSD_before,SSD_after,SSD_after/SSD_before)
%%%%%%%%%%%%%%%

%% Showing the registration result %%%%%
figure;
subplot( 1,3,1);
colormap gray;
imagesc( R );
title('Reference image R')
%     xlabel(['\alpha = ',num2str(alpha)])
axis on;
axis square;

subplot( 1,3,2);
colormap gray;
imagesc( T );
axis on;
axis square;
title('Template image T')
  
subplot( 1,3,3);
colormap gray;
imagesc(Tk); 
title(['Registered Image T_{k}']);
xlabel(['SSDafter = ',num2str(SSD_after),', ReSSD =',...
         num2str(SSD_after/SSD_before)])
axis on;
axis square;
  

