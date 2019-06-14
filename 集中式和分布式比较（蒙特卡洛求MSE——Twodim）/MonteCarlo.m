% Monte_Carlo_Method，进行集中式和分布式估计MSE比较
clc;
clear all;
Len=200;
z0 = 2; % 噪声w的维数
z1 = 1; % 噪声v的维数
to=[-0.15 -0.12 -0.09 -0.05 0 0.05 0.09 0.12 0.15 -0.15 -0.12 -0.09 -0.05 0 0.05 0.09 0.12 0.15];

err_c=zeros(1,Len);
err_d=zeros(19,Len);
err_MHE=zeros(1,Len);
N_Monte=200;

for t1=1:N_Monte
    w=(0.4*rand(z0,Len)-0.2);
    v=(0.4*rand(z1,Len)-0.2);
    [err_1,~,~]=Error_Centralized(w,v,Len,to);
    [err_2,~,~]=Error_Distributed(w,v,Len,to);
    [err_3,~,~]=Error_MHE(w,v,Len,to);
    err_c=err_c+err_1;
    err_d=err_d+err_2;
    err_MHE=err_MHE+err_3;
end
err_c=err_c./N_Monte;
err_d=err_d./N_Monte;
err_MHE=err_MHE./N_Monte;

figure(1)
plot(50:Len,err_c(50:Len),'r-*','LineWidth',2);
hold on
plot(50:Len,err_d(18,50:Len),'b-*','LineWidth',2);
h1=legend({'CEF','DEF'},'FontSize',12)
set(h1,'Orientation','horizon')
xlabel('Time(t)','FontSize',12);
ylabel('MSE','FontSize',12);
axis([50 200 0 0.012]);
set(gcf,'Pos',[300,100,800,250]);

figure(2)
plot(102:Len,err_c(102:Len),'r-*','LineWidth',2);
hold on
plot(102:Len,err_MHE(102:Len),'g-.','LineWidth',2);
h2=legend({'CEF','MHE'},'FontSize',12)
set(h2,'Orientation','horizon')
xlabel('Time(t)','FontSize',12);
ylabel('MSE','FontSize',12);
axis([102 200 0 0.004]);
set(gcf,'Pos',[300,100,800,250]);

figure(3)
E(19)=plot(50:Len,err_d(19,50:Len),'-*','LineWidth',2); 
hold on
for ss=1:18
    E(ss)=plot(50:Len,err_d(ss,50:Len),'-.','LineWidth',2);  
end
xlabel('Time(t)','FontSize',12);
ylabel('MSE','FontSize',12);
legend([E(19) E(1) E(18)],{'DFE','Local estimator 1','Local estimator 17'},'FontSize',12,'Orientation','horizontal')
% legend([E(18) E(1) E(2) E(3) E(4) E(5) E(6) E(7) E(8) E(9) E(10) E(11) E(12) E(13) E(14) E(15) E(16) E(17)],{'Distributed estimator','Local estimator 1','Local estimator 1','Local estimator 1','Local estimator 1','Local estimator 1','Local estimator 1','Local estimator 1','Local estimator 1','Local estimator 1','Local estimator 1','Local estimator 1','Local estimator 1','Local estimator 1','Local estimator 1','Local estimator 1','Local estimator 1','Local estimator 1'},'Orientation','horizontal')
set(gcf,'Pos',[300,100,800,250]);
axis([50 200 0 0.011]);