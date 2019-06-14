%% MHE
close all
clc
clear

Len=300;
%% 系统参数
A=[0.9021   0.09671   0.04876  0.001642;
   -1.918    0.9021    0.9507   0.04876;
   0.04876  0.001642    0.9508   0.09835;
   0.9507   0.04876   -0.9671    0.9508];
B=[0.1   0   0   0;
     0 0.1   0   0;
     0   0 0.1   0;
     0   0   0 0.1];
Ci =[0 0 1 0];
Di=[0.1 0.1 0.1 0.1];
n = 4; % 状态x的维数
m = size(Ci,1);
L = size(Ci,2)/n; % 二值传感器数目
z0 = 4; % 噪声w的维数
z1 = size(Di,2)/L; % 噪声v的维数
%% 阈值设置
to=[0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8 -0.9];
  %% 噪声设置
w=(0.7*rand(z0,Len)-0.35);
v=(0.4*rand(z1,Len)-0.2);
%% 系统初始化
x(:,1)=[0.618;0;1;0]; % 真实状态初始化
z(1,1:L)=2; % 每个二值传感器的进行比较的量
y(1,1:L)=ones(1,L); % 阈值比较结果初始化
index_switch(1,1:L)=ones(1,L); % 切换时刻index
%% MHE初始化
x_MHE(1:n,1)=[0.618;0;1;0]; % 初始化估计器
err(1,1)=1; % 初始化误差
N=100;
Q=1*eye(4);
R=1;
P=0.00008*eye(4);
for i=2:1:Len
    x((i-1)*n+1:i*n,1)=A*x((i-2)*n+1:(i-1)*n,1)+B*w(:,i);
    for j=1:1:L
        C=Ci(1,(j-1)*n+1:j*n);
        D=Di(1,(j-1)*z1+1:j*z1);
        z(i,j)=C*x((i-1)*n+1:i*n,1)+D*v(:,i);
        y(i,j)=sign(z(i,j)-to(j));
        if(y(i,j)*y(i-1,j)<0)
            index_switch(i,j)=1;
            xi(:,:,i,j)=C'*R*C;
            pi(:,:,i,j)=C'*R*to(j);
        else
            index_switch(i,j)=0;
            xi(:,:,i,j)=zeros(n,n);
            pi(:,:,i,j)=zeros(n,1);
        end
    end
    index=index_switch(i,:);
    term=find(index==1);
    L_index=length(term);
    xi_tall(1:n,1:n,i)=zeros(n,n);
    pi_tall(1:n,1,i)=zeros(n,1);
    for m=1:L
        xi_tall(1:n,1:n,i)=xi_tall(1:n,1:n,i)+xi(1:n,1:n,i,m);
        pi_tall(1:n,1,i)=pi_tall(1:n,1,i)+pi(1:n,1,i,m);
    end
    %% MHE
    if(i<=N+1)
        x_MHE(1:n,i)=x((i-1)*n+1:i*n,1);
    else
        M=kron(diag(repmat(1,1,N),1),-A'*Q)+kron(diag(repmat(1,1,N),-1),-Q*A)+kron(diag(repmat(1,1,N+1)),A'*Q*A);
        MM(1:n,1:n)=P+xi_tall(:,:,i-N+1);
        for j=2:N
            MM(n*(j-1)+1:n*j,n*(j-1)+1:n*j)=Q+xi_tall(:,:,i-(N-j));
        end
        MM(n*N+1:n*N+n,n*N+1:n*N+n)=Q-A'*Q*A;
        M=M+MM;
        
        U(1:n,1)=P*A*x_MHE(1:n,i-N-1)+pi_tall(:,:,i-N+1);
        for j=2:N
            U(n*(j-1)+1:n*j,1)=pi_tall(:,:,i-(N-j));
        end
        U(n*N+1:n*N+n,1)=zeros(n,1);
        x_MHE(1:n,i-N:i)=reshape(pinv(M)*U,n,N+1);
    end
end

x1=x(1:n:n*Len,1);
x1_MHE=x_MHE(1,:);

figure(1)
plot(1:1:Len,x1(1:Len),'-.r','LineWidth',2);
hold on
% plot(1:1:Len,xe11,'-k','LineWidth',1.5);
% plot(1:1:Len,xe21,'-g','LineWidth',1.5);
plot(1:1:Len-1,x1_MHE(1:Len-1),'-b','LineWidth',2);
% plot(1:1:Len,z(:,1),'xc','LineWidth',1.5);
% plot(20:1:Len,z(20:Len,1),'c','LineWidth',1.5);
legend({'x','Fusion Estimation of x'},'FontSize',12)
xlabel('Time，t');
ylabel('State');
% axis([50 200 0.1 1.2]);
set(gcf,'Pos',[300,100,700,400]);

x2=x(1:n:n*Len,1);
x2_MHE=x_MHE(1,:);
figure(2)
plot(1:1:Len,x2(1:Len),'-.r','LineWidth',2);
hold on
% plot(1:1:Len,xe11,'-k','LineWidth',1.5);
% plot(1:1:Len,xe21,'-g','LineWidth',1.5);
plot(1:1:Len-1,x2_MHE(1:Len-1),'-b','LineWidth',2);
% plot(1:1:Len,z(:,1),'xc','LineWidth',1.5);
% plot(20:1:Len,z(20:Len,1),'c','LineWidth',1.5);
legend({'x','Fusion Estimation of x'},'FontSize',12)
xlabel('Time，t');
ylabel('State');
% axis([50 200 0.1 1.2]);
set(gcf,'Pos',[300,100,700,400]);

% figure(5)
% plot(err(1,:),'b','LineWidth',1.5);
% hold on
% plot(err(2,:),'g','LineWidth',1.5);
% plot(err(3,:),'k','LineWidth',1.5);
% % plot(err(4,:),'k','LineWidth',1.5);


