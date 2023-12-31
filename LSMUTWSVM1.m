%leastsquarse+ multi-tasks(3-task)+UNIVERSUM TWSVM1:

function [u0,u1,u2,u3]=LSMUTWSVM1(A1,A2,A3,B1,B2,B3,U1,U2,U3,epsilon,T,mu1,C,C_u)

A=[A1; A2; A3];
B=[B1; B2; B3];
U=[U1; U2; U3];

[m1,n]=size(A);
[m11,n]=size(A1);
[m12,n]=size(A2);
[m13,n]=size(A3);
% [m14,n]=size(A4);
[m2,~]=size(B);
[m21,~]=size(B1);
[m22,~]=size(B2);
[m23,~]=size(B3);
% [m24,~]=size(B4);
[u,~]=size(U);
[uu1,~]=size(U1);
[uu2,~]=size(U2);
[uu3,~]=size(U3);
% [uu4,~]=size(U4);
e1=ones(m1,1);
e2=ones(m2,1);
e_u=ones(u,1);
%%%%%%%
A=[A e1];
B=[B e2];
U=[U e_u];
A1=[A1 ones(m11,1)];
A2=[A2 ones(m12,1)];
A3=[A3 ones(m13,1)];
% A4=[A4 ones(m14,1)];
B1=[B1 ones(m21,1)];
B2=[B2 ones(m22,1)];
B3=[B3 ones(m23,1)];
% B4=[B4 ones(m24,1)];
U1=[U1 ones(uu1,1)];
U2=[U2 ones(uu2,1)];
U3=[U3 ones(uu3,1)];
% U4=[U4 ones(uu4,1)];
%%%%%%%%%%%%%%
Q1=B*(inv(A'*A+10^-4*eye(7)))*B';
Q2=B*(inv(A'*A+10^-4*eye(7)))*U';
P11=B1*(inv(A1'*A1+10^-4*eye(7)))*B1';
P12=B2*(inv(A2'*A2+10^-4*eye(7)))*B2';
P13=B3*(inv(A3'*A3+10^-4*eye(7)))*B3';
% P14=B4*(inv(A4'*A4))*B4';
P21=B1*(inv(A1'*A1+10^-4*eye(7)))*U1';
P22=B2*(inv(A2'*A2+10^-4*eye(7)))*U2';
P23=B3*(inv(A3'*A3+10^-4*eye(7)))*U3';
% P24=B4*(inv(A4'*A4))*U4';
P1=blkdiag(P11,P12,P13);
P2=blkdiag(P21,P22,P23);
S1=U*(inv(A'*A+10^-4*eye(7)))*B';
S2=U*(inv(A'*A+10^-4*eye(7)))*U';
R11=U1*(inv(A1'*A1+10^-4*eye(7)))*B1';
R12=U2*(inv(A2'*A2+10^-4*eye(7)))*B2';
R13=U3*(inv(A3'*A3+10^-4*eye(7)))*B3';
% R14=U4*(inv(A4'*A4))*B4';
R21=U1*(inv(A1'*A1+10^-4*eye(7)))*U1';
R22=U2*(inv(A2'*A2+10^-4*eye(7)))*U2';
R23=U3*(inv(A3'*A3+10^-4*eye(7)))*U3';
% R24=U4*(inv(A4'*A4))*U4';
R1=blkdiag(R11,R12,R13);
R2=blkdiag(R21,R22,R23);
I1=eye(m2);
I2=eye(u);
H1=[Q1 -Q2;-S1 S2]+(T/mu1)*[P1 -P2;-R1 R2]+[(1/C)*I1 zeros(m2,u);zeros(u,m2) (1/C_u)*I2];
H2=[e2;(-1+epsilon)*e_u];
x=inv(H1)*H2;
alpha=x(1:m2);
beta=x(m2+1:m2+u);
u0=-inv(A'*A+10^-4*eye(7))*(B'*alpha-U'*beta);
u1=-(T/mu1)*inv(A1'*A1+10^-4*eye(7))*(B1'*alpha(1:m21)-U1'*beta(1:uu1));
u2=-(T/mu1)*inv(A2'*A2+10^-4*eye(7))*(B2'*alpha(m21+1:m21+m22)-U2'*beta(uu1+1:uu1+uu2));
u3=-(T/mu1)*inv(A3'*A3+10^-4*eye(7))*(B3'*alpha(m21+m22+1:m21+m22+m23)-U3'*beta(uu1+uu2+1:uu1+uu2+uu3));
% u4=-(T/mu1)*inv(A4'*A4)*(B4'*alpha(m21+m22+m23+1:m21+m22+m23+m24)-U2'*beta(uu1+uu2+uu3+1:uu1+uu2+uu3+uu4));


