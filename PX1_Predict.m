%% 本程序实现Assumed Density Filtering
% Inputs：
% Outputs：
% Main Steps:
%        1-Prediction
%        2-Update
%        3-Collapse
clc;clear;close all
Q=0.001;% 过程噪声,位移量
R=0.01;% 测量噪声,位移量
C=[1 0];vs=0.8;%行人的速度
% 坐标网格化
nx=0.1;ny=0.1;a=-2:nx:2;b=-2:ny:2;
[X,Y]=meshgrid(a,b);%坐标网格化,X-位移；Y-速度
%% Step1-Prediction， 
X0_Mu0=[0 vs];X0_Sigma0=[0.01 1e-4;1e-4 0.04];% X0的先验分布参数
S=2;A=TransMatrix(S);% S决定A
PdfData=[];
for m=1:size(X,2) % 外两层循环改变X1的取值    
    for n=1:size(Y,1)
        x1=[X(1,m) Y(n,1)];   
        Temp3=0;% 针对每一个X1，积分结果保存在Temp3中
        for i=1:size(X,2)% 内两层循环改变X0的取值，对X0积分
            for j=1:size(Y,1)
                xr=[X(1,i) Y(j,1)];%积分变量
                Temp1=mvnpdf(xr,X0_Mu0,X0_Sigma0);
                X1_Mu=(A*xr')';X1_Sigma=X0_Sigma0;
                Temp2=mvnpdf(x1,X1_Mu,X1_Sigma);
                Temp3=Temp3+Temp1*Temp2;
            end
        end
        PdfData=[PdfData;m n Temp3];%PdfData的前两列表示点的坐标index，第三列为概率密度
        
    end
end
% 当概率密度函数的积分不为1，需要归一化操作
Alpha=sum(PdfData(:,3))*(nx*ny);
PdfData(:,3)=PdfData(:,3)/Alpha;%概率密度函数归一化

figure
subplot(2,1,1)
contour(X,Y,reshape(PdfData(:,3),size(X)))
grid on
xlabel('Lateral Position [m]')
ylabel('Velocity [m/s]')
set(gca,'ylim',[0 3])
set(gca,'ytick',0:0.5:3)
title('True Distribution-P(X1|S1)')


[Mux,Sigx,Muy,Sigy,Covxy]=EstimateBinGaussianPara(PdfData,X,Y);%估计高斯分布的参数
Mu_a=[Mux Muy];%均值
Sigma_a=[Sigx Covxy;Covxy Sigy];%方差
subplot(2,1,2)
contour(X,Y,reshape(mvnpdf([X(:) Y(:)],Mu_a,Sigma_a),size(X)))
grid on
xlabel('Lateral Position [m]')
ylabel('Velocity [m/s]')
set(gca,'ylim',[0 3])
set(gca,'ytick',0:0.5:3)
title('Approximated Distribution-P(X1|S1)')

save Mu_P Mu_a
save Sigma_P Sigma_a
