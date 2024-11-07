function [E,h]=check_B02(W,E0,h,Edot,threshold,epsil)
hcorr=0;hnuovo=0;n=length(E0);norm(Edot,'fro')
Eh = E0 + h*Edot;Eh=Eh/norm(Eh,'fro');
Bt = W+epsil*Eh;P0=Bt<0;P0=P0.*ones(n);

while hcorr==0
disp('correzione h')
Eh = E0 + h*Edot;Eh=Eh/norm(Eh,'fro');E0=h*E0/norm(E0,'fro');


    %beta=norm(Eh,'fro');
    %E = Eh/beta;
    %   beta changes E which implies that some entry of B may become negative
    Bt = W+epsil*Eh;conta=0;
    
conta=1;PIn=zeros(n);P1=W>0;P00=zeros(n);
if min(min(Bt))==0
    hcorr=1;
end
while min(min(Bt))<0 || hcorr==0%&& abs(norm(Eh,'fro')-1)>threshold
P1=W>0;PIn=zeros(n);P1=P1-P0;
min(min(Bt))
if min(min(Bt))>-threshold
    
    E=Eh;
    return
end
    for i=1:n
for j=1:i-1
if P0(i,j)==1
    
theta=-(W(i,j)+epsil*E0(i,j))/(epsil*h*Edot(i,j));

if theta<=-1 
hnuovo=1;conta=conta+1;
end
conta=conta+1;
Eh(i,j)=E0(i,j)+theta*h*Edot(i,j);
Eh(j,i)=E0(i,j)+theta*h*Edot(i,j);
P0(i,j)=1;P0(j,i)=1;
P1(i,j)=0;P1(j,i)=0;
end
end
end

if norm(P0.*Eh,'fro')>1-threshold

rho2=sqrt(1-norm(P0.*E0,'fro')^2)/norm(P1.*Eh,'fro');
else

rho2=sqrt(1-norm(P0.*Eh,'fro')^2)/norm(P1.*Eh,'fro');
end
Eh=P0.*Eh+rho2*P1.*Eh;
Bt=W+epsil*Eh;
for i=1:n
for j=1:i-1
if Bt(i,j)<0
Eh(i,j)=-W(i,j)/epsil;Eh(j,i)=Eh(i,j);
PIn(i,j)=1;PIn(j,i)=1;
P1(i,j)=0;P1(j,i)=0;%P0(i,j)=1;P0(j,i)=1;
end
end
end
P00=(PIn+P0)>0;
if norm(P00.*Eh,'fro')>1-threshold

rho2=sqrt(1-norm(P00.*E0,'fro')^2)/norm(P1.*Eh,'fro');
else
rho2=sqrt(1-norm(P00.*Eh,'fro')^2)/norm(P1.*Eh,'fro');
end
Eh=P00.*Eh+rho2*P1.*Eh;
%Eh=Eh/norm(Eh,'fro');
Bt=W+epsil*Eh;

end

if hnuovo==0
hcorr=1;

else
h=2*h
hnuovo=[];hnuovo=0;

end

%hcorr=1;
end
disp('vediamo')
norm(Eh,'fro')
abs(norm(Eh,'fro')-1);
norm(Eh,'fro');
E=Eh;%E=E/norm(E,'fro');