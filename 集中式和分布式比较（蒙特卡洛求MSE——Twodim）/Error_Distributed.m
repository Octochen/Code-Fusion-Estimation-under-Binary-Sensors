function [err,x_fusion,x] = Error_Distributed(w,v,Len,to)
% 求取估计误差
% 分布式
% 一维时变系统
% 有界噪声方法（特别注意噪声的设置）

Ci=[1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];
% Ci=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
Di = [0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2];
n = 2; % 状态x的维数
m = size(Ci,1);
L = size(Ci,2)/n; % 二值传感器数目

z0 = 2; % 噪声w的维数
z1 = size(Di,2)/L; % 噪声v的维数

x(1:n,1)=[2;2]; % 真实状态初始化
z(1,1:L)=2; % 每个二值传感器的进行比较的量
xe(:,1:L)=2*ones(2,L); % 初始化所有估计器
y(1,1:L)=ones(1,L); % 阈值比较结果初始化
index_switch(1,1:L)=ones(1,L); % 切换时刻index
dt(1,1:L)=zeros(1,L); % 每个时刻补偿步长
Weight=eye(n)/L;
for k=2:1:L
    Weight=[Weight eye(n)/L];
end
W(1:n,1:n*L)=Weight; % 初始化融合权重
x_fusion(1:n,1)=[2;2]; % 初始化融合估计
err=zeros(L+1,Len);
err(1:L+1,1)=ones(L+1,1); % 初始化误差

for i=2:1:Len
    x((i-1)*n+1:i*n,1)=A(i-1)*x((i-2)*n+1:(i-1)*n,1)+B(i-1)*w(:,i);
    for j=1:1:L    
        C=Ci(1,(j-1)*n+1:j*n);
        D=Di(1,(j-1)*z1+1:j*z1); 
        z(i,j)=C*x((i-1)*n+1:i*n,1)+D*v(:,i);
        y(i,j)=sign(z(i,j)-to(j));
        K_switch(:,j)=Bounded_Noises_Local_Gain(A(i-1),B(i-1),C,D);
        if(y(i,j)*y(i-1,j)<0)
            index_switch(i,j)=1;
            dt(i,j)=0;
            xe((i-1)*n+1:i*n,j)=A(i-1)*xe((i-2)*n+1:(i-1)*n,j)+K_switch(:,j)*[to(j)-0.5*C*A(i-1)*xe((i-2)*n+1:(i-1)*n,j)-0.5*C*xe((i-2)*n+1:(i-1)*n,j)];
        else
            index_switch(i,j)=0;
            xe((i-1)*n+1:i*n,j)=A(i-1)*xe((i-2)*n+1:(i-1)*n,j);
        end
        err(j,i)=[x((i-1)*n+1:i*n,1)-xe((i-1)*n+1:i*n,j)]'*[x((i-1)*n+1:i*n,1)-xe((i-1)*n+1:i*n,j)];
    end
    
    for k=1:1:L
        x_total((k-1)*n+1:k*n,:)=xe((i-1)*n+1:i*n,k);
    end
    
    index=index_switch(i,:);
    if(sum(index)==0)
        W((i-1)*n+1:i*n,:)=W((i-2)*n+1:(i-1)*n,:);
    elseif(sum(index)==1)
        W((i-1)*n+1:i*n,:)=zeros(n,n*L);
        ind_k=find(index==1);
        W((i-1)*n+1:i*n,(ind_k-1)*n+1:ind_k*n)=eye(n);
    elseif(sum(index)>1)
        term=find(index==1);
        L_term=length(term);
        
        A_F=[];
        K_M=[];
        B_M=[];
        C_up=[];
        C_down=[];
        D_up=[];
        D_down=[];
        for k=1:1:L_term
            C=Ci(1,(term(k)-1)*n+1:term(k)*n);
            D=Di(1,(term(k)-1)*m+1:term(k)*m);  
            K_M=mdiag(K_M,K_switch(:,term(k)));
            B_M=[B_M;B(i-1)];
            C_up=[C_up;0.5*C];
            C_down=[C_down;0.5*C];
            D_up=mdiag(D_up,0.5*D);
            D_down=mdiag(D_down,0.5*D);
            
            A_F_term=A(i-1)-0.5*K_switch(:,term(k))*C*A(i-1)-0.5*K_switch(:,term(k))*C;
            A_F=mdiag(A_F,A_F_term);
        end
        B_F_ba=[zeros(n*L_term,n) B_M-K_M*C_up*B(i-1) -K_M*D_up -K_M*D_down];
        
        setlmis([])
            eps=lmivar(1,[1,1]);
            P=lmivar(1,[n*L_term,1]);  
            Theta=lmivar(1,[n+z0+2*L_term*z1,1]);
            Gamma=lmivar(2,[n*L_term,n+z0+2*L_term*z1]);
            [W_term1,nW_term1,sW_term1]=lmivar(2,[n,n*(L_term-1)]);
            [W_term2,nW_term2,sW_term2]=lmivar(3,[sW_term1 zeros(n)]);
    
            I_term=[];
            for I_term_m=1:1:L_term
                I_term=[I_term;-eye(n)];
            end
            I_total=[zeros(n*L_term,n*(L_term-1)) I_term];
    
            %1st LMI
            lmiterm([1 1 1 eps],-1,1)
            lmiterm([1 1 2 -W_term2],K_M',1)
            
            lmiterm([1 2 2 0],-1);
            lmiterm([1 2 3 W_term2],1,A_F); 
            lmiterm([1 2 3 0],[zeros(n,n*(L_term-1)) eye(n)]*A_F);
            lmiterm([1 2 3 W_term2],1,I_total*A_F);
                                  
            lmiterm([1 2 4 W_term2],1,B_F_ba);
            lmiterm([1 2 4 0],[zeros(n,n*(L_term-1)) eye(n)]*B_F_ba);
            lmiterm([1 2 4 W_term2],1,I_total*B_F_ba);
            
            lmiterm([1 3 3 P],-1,1);
            lmiterm([1 3 4 Gamma],-1,1);    
            
            lmiterm([1 4 4 Theta],-1,1);
            lmiterm([1 4 5 eps],1,[C_down-C_up*A(i-1) -C_up*B(i-1) -D_up -D_down]');
            
            lmiterm([1 5 5 eps],-1,1);    
            
            
            lmisys=getlmis;
        %LMI 求解
        n_LMI=decnbr(lmisys);
        c_LMI=zeros(1,n_LMI);
        for j=1:n_LMI
            [P_j,Theta_j]=defcx(lmisys,j,P,Theta);
            c_LMI(j)=trace(P_j)+trace(Theta_j);
        end  
        [copt xopt]=mincx(lmisys,c_LMI); % 求可行解
        W_opt_term2=dec2mat(lmisys,xopt,W_term2)  
        
        W_opt=zeros(n,n*L);
        for k=1:1:L_term-1
            W_opt(:,(term(k)-1)*n+1:term(k)*n)=W_opt_term2(:,(k-1)*n+1:k*n);
        end
        W_opt(:,(term(L_term)-1)*n+1:term(L_term)*n)=eye(n)+W_opt_term2*I_term;        
        W((i-1)*n+1:i*n,:)=W_opt;
    end
    
    x_fusion(1:n,i)=W(n*(i-1)+1:n*i,:)*x_total;
    err(L+1,i)=[x((i-1)*n+1:i*n,1)-x_fusion(1:n,i)]'*[x((i-1)*n+1:i*n,1)-x_fusion(1:n,i)];
end

end

