function P=pseudoinv(L,X,lam,ind)
% L = n dim Laplacian matrix
% X = (n x (at least ind)) dim matrix of the eigenvectors from 1st to ind-th
% lam = eigenvalue
% ind = numbers of eigenvectors


dim=length(L(:,1));
X1=X;
X1(:,ind)=X(:,1);
X1(:,1)=X(:,ind);
dd=X1'*(L-lam*eye(dim))*X1;
d=1./diag(dd);
D=diag([0;d(2:end)]);
P=X1*D*X1';
P=(P+P')/2;
