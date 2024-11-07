function [v,mean_vP]=constr_matrix_m(vett_mem,x2)
% vett_mem is a vector of entries 0,1,-1 (in accordance with membership constraint)
ok1=sum((vett_mem.*x2)>0);
ok2=sum((-vett_mem.*x2)>0);
if ok2>ok1
    vett_mem=-vett_mem;
end
dim=length(vett_mem);
oneP=vett_mem>0;oneP=ones(dim,1).*oneP;
oneM=vett_mem<0;oneM=ones(dim,1).*oneM;

mean_vP=mean(x2.*oneP')*dim/sum(oneP);
mean_vM=mean(x2.*oneM')*dim/sum(oneM);
v=0;
for i=1:dim
    vectE=zeros(dim,1);vectE(i)=1;
    if vett_mem(i)>0
        v=v-(x2(i)-mean_vP)*(vectE-oneP./sum(oneP));
        Xmean(i)=mean_vP;
    elseif vett_mem(i)<0
        v=v-(x2(i)-mean_vM)*(vectE-oneM./sum(oneM));
        Xmean(i)=mean_vM;
    end
end