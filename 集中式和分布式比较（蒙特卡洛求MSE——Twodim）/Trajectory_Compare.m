% 集中式和分布式融合迹比较
clc;
clear
close all;
Len=200;
z0 = 2; % 噪声w的维数
z1 = 1; % 噪声v的维数
w=(0.4*rand(z0,Len)-0.2);
v=(0.4*rand(z1,Len)-0.2);
to=[-0.15 -0.12 -0.09 -0.05 0 0.05 0.09 0.12 0.15 -0.15 -0.12 -0.09 -0.05 0 0.05 0.09 0.12 0.15];
% to=[-0.4 -0.36 -0.32 -0.27 -0.22 -0.17 -0.12 -0.07 -0.02 0.02 0.7 0.12 0.17 0.22 0.27 0.32 0.36 0.4];

[err_c,fusion_c,x_c]=Error_Centralized(w,v,Len,to);
[err_d,fusion_d,x_d]=Error_Distributed(w,v,Len,to);
[err_MHE,fusion_MHE,x_MHE]=Error_MHE(w,v,Len,to);

n=2;
x_d1=x_d(1:n:n*Len,1);
x_d2=x_d(2:n:n*Len,1);
x_fusion_d1=fusion_d(1,50:Len);
x_fusion_d2=fusion_d(2,50:Len);
x_fusion_c1=fusion_c(1,50:Len);
x_fusion_c2=fusion_c(2,50:Len);
x_fusion_MHE1=fusion_MHE(1,50:Len);
x_fusion_MHE2=fusion_MHE(2,50:Len);

figure(1)
plot(50:Len,x_d1(50:Len),'k','LineWidth',2);
hold on
plot(50:Len,x_fusion_c1,'r-.','LineWidth',2);
% plot(50:Len,x_fusion_MHE1,'g-.','LineWidth',2);
legend({'x','x_{CFE}'},'FontSize',12)
xlabel('Time(t)','FontSize',12);
ylabel('State 1','FontSize',12);
set(gcf,'Pos',[300,100,800,250]);

figure(2)
plot(50:Len,x_d2(50:Len),'k','LineWidth',2);
hold on
plot(50:Len,x_fusion_c2,'r-.','LineWidth',2);
% plot(50:Len,x_fusion_MHE2,'g-.','LineWidth',2);
legend({'x','x_{CFE}'},'FontSize',12)
xlabel('Time(t)','FontSize',12);
ylabel('State 2','FontSize',12);
set(gcf,'Pos',[300,100,800,250]);

figure(3)
plot(50:Len,x_d1(50:Len),'k','LineWidth',2);
hold on
plot(50:Len,x_fusion_d1,'b-.','LineWidth',2);
% plot(50:Len,x_fusion_MHE1,'g-.','LineWidth',2);
legend({'x','x_{DFE}'},'FontSize',12)
xlabel('Time(t)','FontSize',12);
ylabel('State 1','FontSize',12);
set(gcf,'Pos',[300,100,800,250]);

figure(4)
plot(50:Len,x_d2(50:Len),'k','LineWidth',2);
hold on
plot(50:Len,x_fusion_d2,'b-.','LineWidth',2);
% plot(50:Len,x_fusion_MHE2,'g-.','LineWidth',2);
legend({'x','x_{DFE}'},'FontSize',12)
xlabel('Time(t)','FontSize',12);
ylabel('State 2','FontSize',12);
set(gcf,'Pos',[300,100,800,250]);