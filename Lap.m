function L=Lap(M)
% compute the laplacian matrix L by the adjacency matrix M

n=length(M);
One=ones(n,1);
L=diag(M*One)-M;

return
end

