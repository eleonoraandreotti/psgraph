clc
close all
clear

h0=1e-5;
T=100;
Wtype=1;

%----karate example
W=karate();
d=length(W);
soglia=17;x_fix=zeros(d,1);
%----second exemple
%     riga=zeros(d,1);riga(2:3)=1*ones(2,1);
%     W=toeplitz(riga);W(1,3)=0;W(3,1)=0;
%     d=8;
%---------

%  3 dim vector of coefficients:
% 1st entry= for cutting (usually =1)
% 2nd entry= cardinality constraint
% 3rdentry= membership constraint
Cost=[1,10,0]; 


kmax=20;
eps0=3;
eps_lb=1;eps_ub=10;
tol=1e-6;
[epsil,E,delta] = compclos( W,h0,T,kmax,tol,Cost,soglia,x_fix,eps0,eps_lb,eps_ub);
B=W+epsil*E;
[V1,V2]=eig(Lap(B));v1=V1(:,1);v2=V1(:,2);
