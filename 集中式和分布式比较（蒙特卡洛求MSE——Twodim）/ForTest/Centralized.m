% 集中式
% 2维时变系统
% 有界噪声方法（特别注意噪声的设置）

close all
clc
clear

Len=200;

Ci=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
Di = [0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2];
% to=[0.10 0.18 0.26 0.34 0.42 0.5 0.58 0.66 0.74 0.82 0.9 0.98]; % 所有二值传感器的阈值
to=[-0.4 -0.35 -0.3 -0.25 -0.2 -0.15 -0.1 -0.05 0 0.05 0.10 0.15 0.2 0.25 0.3 0.35 0.4];
n = 2; % 状态x的维数
m = size(Ci,1);
L = size(Ci,2)/n; % 二值传感器数目

z0 = 2; % 噪声w的维数
z1 = size(Di,2)/L; % 噪声v的维数

w=(0.7*rand(z0,Len)-0.35);
v=(0.4*rand(z1,Len)-0.2);

x(1:n,1)=[2;2]; % 真实状态初始化
z(1,1:L)=2; % 每个二值传感器的进行比较的量
y(1,1:L)=ones(1,L); % 阈值比较结果初始化
index_switch(1,1:L)=ones(1,L); % 切换时刻index
dt(1,1:L)=zeros(1,L); % 每个时刻补偿步长
x_fusion(1:n,1)=[2;2]; % 初始化集中式估计
err(1,1)=1; % 初始化误差

for i=2:1:Len
    x((i-1)*n+1:i*n,1)=A(i-1)*x((i-2)*n+1:(i-1)*n,1)+B(i-1)*w(:,i);
    for j=1:1:L
        C=Ci(1,(j-1)*n+1:j*n);
        D=Di(1,(j-1)*z1+1:j*z1);
        z(i,j)=C*x((i-1)*n+1:i*n,1)+D*v(:,i);
        y(i,j)=sign(z(i,j)-to(j));
        if(y(i,j)*y(i-1,j)<0)
            index_switch(i,j)=1;
            dt(i,j)=0;
        else
            index_switch(i,j)=0;
        end
    end
    
    index=index_switch(i,:);
    if(sum(index)==0)
        x_fusion(1:n,i)=A(i-1)*x_fusion(1:n,i-1);
    else
        term=find(index==1);
        L_index=length(term);
   
        C_up=[];
        C_down=[];
        D_up=[];
        D_down=[];
        tou=[];
        for k=1:1:L_index
            C=Ci(1,(term(k)-1)*n+1:term(k)*n);
            D=Di(1,(term(k)-1)*m+1:term(k)*m);  
            tou=[tou;to(term(k))];    
            C_up=[C_up;0.5*C];
            C_down=[C_down;0.5*C];
            D_up=mdiag(D_up,0.5*D);
            D_down=mdiag(D_down,0.5*D);
        end
        
        setlmis([])
            P=lmivar(1,[n,1]);  
            eps=lmivar(1,[1,1]);
            [K,nK,sK]=lmivar(2,[n,L_index]);
            Theta=lmivar(1,[n+z0+2*z1*L_index,1]);
            Gamma=lmivar(2,[n,n+z0+2*z1*L_index]);
            %1st LMI
            lmiterm([1 1 1 eps],-1,1);
            lmiterm([1 1 2 -K],1,1);  
        
            lmiterm([1 2 2 0],-1);
            lmiterm([1 2 3 0],A(i-1));
            lmiterm([1 2 3 K],-1,C_up*A(i-1)+C_down);
            lmiterm([1 2 4 0],[zeros(n) B(i-1) zeros(n,z1*L_index) zeros(n,z1*L_index)]);
            lmiterm([1 2 4 K],1,[zeros(L_index,n) -C_up*B(i-1) -D_up -D_down]);
    
            lmiterm([1 3 3 P],-1,1);
            lmiterm([1 3 4 Gamma],-1,1);
    
            lmiterm([1 4 4 Theta],-1,1);
            lmiterm([1 4 5 eps],1,[C_down-C_up*A(i-1) -C_up*B(i-1) -D_up D_down]');
    
            lmiterm([1 5 5 eps],-1,1);    

        lmisys=getlmis;
        %LMI 求解
        n_LMI=decnbr(lmisys);
        c_LMI=zeros(1,n_LMI);
        for j=1:n_LMI
            [Theta_j,P_j]=defcx(lmisys,j,Theta,P);
            c_LMI(j)=trace(Theta_j)+trace(P_j);
        %     c_LMI(j)=trace(Theta_j);
        end
        [copt xopt]=mincx(lmisys,c_LMI); % 求可行解
        K_opt=dec2mat(lmisys,xopt,K) % 提取解矩阵
        
        x_fusion(1:n,i)=A(i-1)*x_fusion(1:n,i-1)+K_opt*[tou-C_up*A(i-1)*x_fusion(1:n,i-1)-C_down*x_fusion(1:n,i-1)];
    end
%     err(L+1,i)=[x((i-1)*n+1:i*n,1)-x_fusion((i-1)*n+1:i*n,1)]'*[x((i-1)*n+1:i*n,1)-x_fusion((i-1)*n+1:i*n,1)];
end

x1=x(1:n:n*Len,1);
x1_fusion=x_fusion(1,:);

figure(1)
plot(50:1:Len,x1(50:Len),'-.r','LineWidth',2);
hold on
% plot(1:1:Len,xe11,'-k','LineWidth',1.5);
% plot(1:1:Len,xe21,'-g','LineWidth',1.5);
plot(50:1:Len,x1_fusion(50:Len),'-*b','LineWidth',2);
% plot(1:1:Len,z(:,1),'xc','LineWidth',1.5);
% plot(20:1:Len,z(20:Len,1),'c','LineWidth',1.5);
legend({'x','Fusion Estimation of x'},'FontSize',12)
xlabel('Time，t');
ylabel('State');
% axis([50 200 0.1 1.2]);
set(gcf,'Pos',[300,100,700,400]);

x2=x(1:n:n*Len,1);
x2_fusion=x_fusion(1,:);
figure(2)
plot(50:1:Len,x2(50:Len),'-.r','LineWidth',2);
hold on
% plot(1:1:Len,xe11,'-k','LineWidth',1.5);
% plot(1:1:Len,xe21,'-g','LineWidth',1.5);
plot(50:1:Len,x2_fusion(50:Len),'-*b','LineWidth',2);
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


