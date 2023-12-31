%leastsquarse+ multi-tasks(3-task)+UNIVERSUM TWSVM2:

function [v0,v1,v2,v3]=LSMUTWSVM2(A1,A2,A3,B1,B2,B3,U1,U2,U3,epsilon,T,mu2,C,C_u)

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
%%%%%%
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
%%%%%%%%%%%

Q1=A*(inv(B'*B+10^-4*eye(7)))*A';
Q2=A*(inv(B'*B+10^-4*eye(7)))*U';
P11=A1*(inv(B1'*B1+10^-4*eye(7)))*A1';
P12=A2*(inv(B2'*B2+10^-4*eye(7)))*A2';
P13=A3*(inv(B3'*B3+10^-4*eye(7)))*A3';
% P14=A4*(inv(B4'*B4))*A4';
P21=A1*(inv(B1'*B1+10^-4*eye(7)))*U1';
P22=A2*(inv(B2'*B2+10^-4*eye(7)))*U2';
P23=A3*(inv(B3'*B3+10^-4*eye(7)))*U3';
% P24=A4*(inv(B4'*B4))*U4';
P1=blkdiag(P11,P12,P13);
P2=blkdiag(P21,P22,P23);
S1=U*(inv(B'*B+10^-4*eye(7)))*A';
S2=U*(inv(B'*B+10^-4*eye(7)))*U';
R11=U1*(inv(B1'*B1+10^-4*eye(7)))*A1';
R12=U2*(inv(B2'*B2+10^-4*eye(7)))*A2';
R13=U3*(inv(B3'*B3+10^-4*eye(7)))*A3';
% R14=U4*(inv(B4'*B4))*A4';
R21=U1*(inv(B1'*B1+10^-4*eye(7)))*U1';
R22=U2*(inv(B2'*B2+10^-4*eye(7)))*U2';
R23=U3*(inv(B3'*B3+10^-4*eye(7)))*U3';
% R24=U4*(inv(B4'*B4))*U4';
R1=blkdiag(R11,R12,R13);
R2=blkdiag(R21,R22,R23);
I1=eye(m1);
I2=eye(u);
H1=[Q1 -Q2;S1 -S2]+(T/mu2)*[P1 -P2;R1 -R2]+[(1/C)*I1 zeros(m1,u);zeros(u,m1) (-1/C_u)*I2];
H2=[e1;(1-epsilon)*e_u];
x=inv(H1)*H2;
alpha1=x(1:m1);
beta1=x(m1+1:m1+u);
v0=inv(B'*B+10^-4*eye(7))*(A'*alpha1-U'*beta1);
v1=(T/mu2)*inv(B1'*B1+10^-4*eye(7))*(A1'*alpha1(1:m11)-U1'*beta1(1:uu1));
v2=(T/mu2)*inv(B2'*B2+10^-4*eye(7))*(A2'*alpha1(m11+1:m11+m12)-U2'*beta1(uu1+1:uu1+uu2));
v3=(T/mu2)*inv(B3'*B3+10^-4*eye(7))*(A3'*alpha1(m11+m12+1:m11+m12+m13)-U3'*beta1(uu1+uu2+1:uu1+uu2+uu3));
% v4=(T/mu2)*inv(B4'*B4)*(A4'*alpha1(m11+m12+m13+1:m11+m12+m13+m14)-U2'*beta1(uu1+uu2+uu3+1:uu1+uu2+uu3+uu4));


