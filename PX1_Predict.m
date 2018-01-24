%% ������ʵ��Assumed Density Filtering
% Inputs��
% Outputs��
% Main Steps:
%        1-Prediction
%        2-Update
%        3-Collapse
clc;clear;close all
Q=0.001;% ��������,λ����
R=0.01;% ��������,λ����
C=[1 0];vs=0.8;%���˵��ٶ�
% ��������
nx=0.1;ny=0.1;a=-2:nx:2;b=-2:ny:2;
[X,Y]=meshgrid(a,b);%��������,X-λ�ƣ�Y-�ٶ�
%% Step1-Prediction�� 
X0_Mu0=[0 vs];X0_Sigma0=[0.01 1e-4;1e-4 0.04];% X0������ֲ�����
S=2;A=TransMatrix(S);% S����A
PdfData=[];
for m=1:size(X,2) % ������ѭ���ı�X1��ȡֵ    
    for n=1:size(Y,1)
        x1=[X(1,m) Y(n,1)];   
        Temp3=0;% ���ÿһ��X1�����ֽ��������Temp3��
        for i=1:size(X,2)% ������ѭ���ı�X0��ȡֵ����X0����
            for j=1:size(Y,1)
                xr=[X(1,i) Y(j,1)];%���ֱ���
                Temp1=mvnpdf(xr,X0_Mu0,X0_Sigma0);
                X1_Mu=(A*xr')';X1_Sigma=X0_Sigma0;
                Temp2=mvnpdf(x1,X1_Mu,X1_Sigma);
                Temp3=Temp3+Temp1*Temp2;
            end
        end
        PdfData=[PdfData;m n Temp3];%PdfData��ǰ���б�ʾ�������index��������Ϊ�����ܶ�
        
    end
end
% �������ܶȺ����Ļ��ֲ�Ϊ1����Ҫ��һ������
Alpha=sum(PdfData(:,3))*(nx*ny);
PdfData(:,3)=PdfData(:,3)/Alpha;%�����ܶȺ�����һ��

figure
subplot(2,1,1)
contour(X,Y,reshape(PdfData(:,3),size(X)))
grid on
xlabel('Lateral Position [m]')
ylabel('Velocity [m/s]')
set(gca,'ylim',[0 3])
set(gca,'ytick',0:0.5:3)
title('True Distribution-P(X1|S1)')


[Mux,Sigx,Muy,Sigy,Covxy]=EstimateBinGaussianPara(PdfData,X,Y);%���Ƹ�˹�ֲ��Ĳ���
Mu_a=[Mux Muy];%��ֵ
Sigma_a=[Sigx Covxy;Covxy Sigy];%����
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
