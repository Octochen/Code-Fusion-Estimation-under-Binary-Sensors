% ���ع�������
% ��ȷ��
function [K_opt] = Bounded_Noises_Local_Gain(A,B,C,D)
n = size(A,1); % ״̬x��ά��
m = size(C,1); % �۲�y_i��ά�� 
z0 = size(B,2); % ����w��ά��
z1 = size(D,2); % ����v_i��ά��
% A=A(t-1)  B=B(t-1)  C=C(t)  D=D(t)

%% ���K
setlmis([])
    P=lmivar(1,[n,1]);  
    eps=lmivar(1,[1,1]);
    [K,nK,sK]=lmivar(2,[n,m]);
    Theta=lmivar(1,[n+z0+2*z1,1]);
    Gamma=lmivar(2,[n,n+z0+2*z1]);
    %1st LMI
    lmiterm([1 1 1 eps],-1,1);
    lmiterm([1 1 2 -K],1,1);
    
    lmiterm([1 2 2 0],-1);
    lmiterm([1 2 3 0],A);
    lmiterm([1 2 3 K],1,-0.5*C-0.5*C*A);
    lmiterm([1 2 4 0],[zeros(n) B zeros(n,z1) zeros(n,z1)]);
    lmiterm([1 2 4 K],1,[zeros(m,n) -0.5*C*B -0.5*D -0.5*D]);
    
    lmiterm([1 3 3 P],-1,1);
    lmiterm([1 3 4 Gamma],-1,1);
    
    lmiterm([1 4 4 Theta],-1,1);
    lmiterm([1 4 5 eps],1,1/2*[C-C*A C*B -D D]');
    
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
end

