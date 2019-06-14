function [err,x_fusion,x] = Error_Centralized(w,v,Len,to)
% ��ȡ�������
% ����ʽ
% һάʱ��ϵͳ
% �н������������ر�ע�����������ã�

Ci=[1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];
% Ci=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
Di = [0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2];
n = 2; % ״̬x��ά��
m = size(Ci,1);
L = size(Ci,2)/n; % ��ֵ��������Ŀ

z0 = 2; % ����w��ά��
z1 = size(Di,2)/L; % ����v��ά��

%% ϵͳ��ʼ��
x(1:n,1)=[2;2]; % ��ʵ״̬��ʼ��
z(1,1:L)=2; % ÿ����ֵ�������Ľ��бȽϵ���
y(1,1:L)=ones(1,L); % ��ֵ�ȽϽ����ʼ��
index_switch(1,1:L)=ones(1,L); % �л�ʱ��index

%% ����ʽ�ںϳ�ʼ��
x_fusion(1:n,1)=[2;2]; % ��ʼ������ʽ����
err(1,1)=1; % ��ʼ�����

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
        %LMI ���
        n_LMI=decnbr(lmisys);
        c_LMI=zeros(1,n_LMI);
        for j=1:n_LMI
            [Theta_j,P_j]=defcx(lmisys,j,Theta,P);
            c_LMI(j)=trace(Theta_j)+trace(P_j);
        %     c_LMI(j)=trace(Theta_j);
        end
        [copt xopt]=mincx(lmisys,c_LMI); % ����н�
        K_opt=dec2mat(lmisys,xopt,K) % ��ȡ�����
        
        x_fusion(1:n,i)=A(i-1)*x_fusion(1:n,i-1)+K_opt*[tou-C_up*A(i-1)*x_fusion(1:n,i-1)-C_down*x_fusion(1:n,i-1)];
    end
    err(1,i)=[x((i-1)*n+1:i*n,1)-x_fusion(1:n,i)]'*[x((i-1)*n+1:i*n,1)-x_fusion(1:n,i)];
end

end

