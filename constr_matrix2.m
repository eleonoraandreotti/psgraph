function [v,mean_x2P,mean_x2M]=constr_matrix2(x2)
dim=length(x2);
oneP=x2>0;oneP=ones(dim,1).*oneP;
oneM=x2<=0;oneM=ones(dim,1).*oneM;
mean_x2P=mean(x2.*oneP)*dim/sum(oneP);
mean_x2M=mean(x2.*oneM)*dim/sum(oneM);
v=0;
for i=1:dim
    vectE=zeros(dim,1);vectE(i)=1;
    if x2(i)>0
        v=v-(x2(i)-mean_x2P)*(vectE-oneP./sum(oneP));
        Xmean(i)=mean_x2P;
    else %if x2(i)<0
        v=v-(x2(i)-mean_x2M)*(vectE-oneM./sum(oneM));
        Xmean(i)=mean_x2M;
    end
end