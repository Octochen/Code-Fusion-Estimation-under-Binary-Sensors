function [err,x_MHE,x] = Error_MHE(w,v,Len,to)
%% MHE

%% ϵͳ����
Ci=[1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];
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
%% MHE��ʼ��
x_MHE(1:n,1)=[2;2]; % ��ʼ��������
err(1,1)=1; % ��ʼ�����
N=100;
Q=1*eye(2);
R=1;
P=0.00008*eye(2);
for i=2:1:Len
    x((i-1)*n+1:i*n,1)=A(i-1)*x((i-2)*n+1:(i-1)*n,1)+B(i-1)*w(:,i);
    for j=1:1:L
        C=Ci(1,(j-1)*n+1:j*n);
        D=Di(1,(j-1)*z1+1:j*z1);
        z(i,j)=C*x((i-1)*n+1:i*n,1)+D*v(j,i);
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
        M=kron(diag(repmat(1,1,N),1),-A(i-1)'*Q)+kron(diag(repmat(1,1,N),-1),-Q*A(i-1))+kron(diag(repmat(1,1,N+1)),A(i-1)'*Q*A(i-1));
        MM(1:n,1:n)=P+xi_tall(:,:,i-N+1);
        for j=2:N
            MM(n*(j-1)+1:n*j,n*(j-1)+1:n*j)=Q+xi_tall(:,:,i-(N-j));
        end
        MM(n*N+1:n*N+n,n*N+1:n*N+n)=Q-A(i-1)'*Q*A(i-1);
        M=M+MM;
        
        U(1:n,1)=P*A(i-1)*x_MHE(1:n,i-N-1)+pi_tall(:,:,i-N+1);
        for j=2:N
            U(n*(j-1)+1:n*j,1)=pi_tall(:,:,i-(N-j));
        end
        U(n*N+1:n*N+n,1)=zeros(n,1);
        x_MHE(1:n,i-N:i)=reshape(pinv(M)*U,n,N+1);
    end
    err(1,i)=[x((i-1)*n+1:i*n,1)-x_MHE(1:n,i)]'*[x((i-1)*n+1:i*n,1)-x_MHE(1:n,i)];
end

end


