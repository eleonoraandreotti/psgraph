function [ epsil,E,Mis] = compclos( W,h0,T,kmax,tol,Cost,soglia,x_fix,eps0,eps_lb,eps_ub)
% INTUP:
% W adjecency matrix 
% h0 time in
% T time out
% kmax max num of iterations
% tol tollerance
% Cost 3-dim vector of constraints: coeff cut, coeff cardinality and coeff
%                                   membership
% soglia: parameter for the cardinality constraints: n/2
% x_fix: vector n dimensional for memebership constraint
% eps_lb,eps_ub: lower and upper bound for epsil

n=length(W);
threshold=1e-8;

Wtype=1;
P=W>0;
epsil=eps0;

    

[Mis,E,ng,Edot]=graph_constr_lin(W,Wtype,0,1,h0,epsil,soglia,x_fix,Cost);
mis=sum(Mis(end,:));
Et=E;
plotta(1,1:5)=[sum(Mis(end,:)),Mis(end,:),epsil];

k=0;
while k<=kmax



  if mis<tol
      eps_ub=min(eps_ub,epsil);
      epsil1=(eps_lb+eps_ub)/2;
      %disp('bisection')
  else
      eps_lb=max(eps_lb,epsil);
      epsil1=epsil-mis/ng;
      %disp('Newton')
  end
  if epsil1<eps_lb || epsil1>eps_ub
      %disp('epsil not in the interval')
      epsil1=(eps_ub+eps_lb)/2;
     
  end
  if (eps_ub-eps_lb)<tol
      return
      [MIS',EPSIL']
  end
      
   
    [E,h_eps]=check_B02(W,E,.000001,Edot,threshold,epsil1);

    [Mis,E,ng,Edot]=graph_constr_lin(W,1,0,10,h0,epsil1,soglia,x_fix,Cost,0,E);  
    mis=sum(Mis(end,:));
    epsil=epsil1;
k=k+1;
MIS(k)=mis;
EPSIL(k)=epsil;
MISTUTTI(k,:)=Mis(end,:);
% for i=1:n
% for j=1:i-1
% 
% if P(i,j)>0
% if B(i,j)<=(THETA*W(i,j))
% CH(i,j)=1;CH(j,i)=1;
% elseif abs(epsil*E(i,j))<=(THETA*W(i,j))
% CH(i,j)=0;CH(j,i)=0;
% else
% CH(i,j)=2;CH(j,i)=2;
% end
% end
% end
% end
% 
% if CH<2
% disp('ci siamo')
% E=-CH.*W/epsil;
% B=W+epsil*E;
% 
% lambda=eig(Lap(B))
% return
% end
 end
[MISTUTTI,EPSIL']
 
