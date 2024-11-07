function [Mis,E0,ng,Edot]=graph_constr_lin(W,Wtype,t0,tf,h0,epsil,soglia,x_fix,Cost,store,Eguess)
% CONSTRAINTS: if they are ==0, they are not considered
% Cost(1)= ( > 0) for the coalescense
% Cost(2)= ( > 0) cardinality constraint
% Cost(3)= ( > 0) membership constraint (by means of x_fix)

% x_fix vector of the membership constraint: its entries are 0 (no
% contsraint), +1 (vertices in a cluster) o -1 (vertices in the other cluster) 
format long
threshold=1e-8; %-threshold è la soglia per la non negatività
Edotmin=sqrt(threshold);
grid=[];
sol=[];
lambda=[];
% Initialization
hmin=1e-8;%20
n=length(W);dim=n;
E=zeros(n,n);
One=ones(n,1);
P=W>threshold;
% pattern
% if nargin<11
%
% else
% P=Eguess>threshold;
% end

% active set if zero
Z=ones(n,n);
P0=zeros(n,n);
% pattern
if (P~=P') | min(min(W))<0
    disp('non symmetric pattern or negative weights present');
    return
end
gamma=2; % per riscalare il passo temporale
reject=0;
h=h0;
ns=1;
facth=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 10
    store=0;
end

Dcost=(Cost>0);


L=Lap(W);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 11
    
    
    [lamsort, X] = compeigvecs(L, Wtype, 0, zeros(n), 3);
    lam2=lamsort(2);x2=X(:,2);
    XX2=sum((X(:,2)-mean(X(:,2))).^2);
    XX1=sum((X(:,1)-mean(X(:,1))).^2);
    if XX1>XX2
        xx2=X(:,1);
    else
        xx2=X(:,2);
    end
    
    
    Mpi=pseudoinv(L,X,lam2,2);
    F=zeros(n);
    if Dcost(1)>0
        
        F2 = x2.^2*One';
        F2 = (F2+F2')/2 - (x2*x2');
        F = Cost(1)*(F2).*P;
        
    end
    if Dcost(2)>0
        
        z=Mpi*constr_matrix2(xx2);
        xc = Cost(2)*z;%(Mpi,x2I,x2);
        Fc = (xx2.*xc)*One'- (xx2*xc');
        Fc = (Fc+Fc')/2 ;
        
        F=F-Fc.*P;
    end
    if Dcost(3)>0
        
        Xm=constr_matrix_m(x_fix,x2');
        xm = Cost(3)*Mpi*Xm;
        
        Fm = (x2.*xm)*One'- (x2*xm');
        Fm = (Fm+Fm')/2 ;
        F=F+Fm.*P;
        
        
    end
    
    
    % projection;
    F =  F - diag(diag(F));
    G = -F;
    E0=G/norm(G,'fro');
    %E0=check_B02(W,E0,.000001,E0,threshold,epsil);

    
else
    E0=Eguess.*P;
    E0=E0/norm(E0,'fro');
end

B=W+epsil*E0;
its=0;



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Check non-negativity of initial matrix (W è non negativa, ma W+epsil*E0 potrebbe non esserlo)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conta=0;
minB=min(min(B));
if minB<-threshold
B
minB
    disp('B è negativa')
    
    return
end
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % The new initial matrix is %
% % B = W+epsil*E0;           %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------
% Compute Graph Laplacian



L = Lap(B);
[lamsort, X] = compeigvecs(L, Wtype, 0, zeros(n), 3);
lam2=lamsort(2);x2=X(:,2);
% if X(:,2)>0
% x2=X(:,1);lam2=lamsort(1);
% elseif X(:,2)>0
% x2=X(:,1);lam2=lamsort(1);
% else
% x2=X(:,2);lam2=lamsort(2);
% end
Mpi=pseudoinv(L,X,lam2,2);%pinv(L-lam2*eye(n));%Mpi=(Mpi+Mpi')/2;
F=zeros(n);
XX2=sum((X(:,2)-mean(X(:,2))).^2);
        XX1=sum((X(:,1)-mean(X(:,1))).^2);
        if XX1>XX2
            xx2=X(:,1);
        else
            xx2=X(:,2);
        end

if Dcost(1)>0
    
    
    F2 = x2.^2*One';
    F2 = (F2+F2')/2 - (x2*x2');
    
    F = Cost(1)*(F2).*P;
    
    Mis(1,1)=Cost(1)*lam2;
else
    Mis(1,1)=Cost(1);
end
if Dcost(2)>0
    
    
    [Xc,meanX2P,meanX2M]=constr_matrix2(xx2);
    xc = Cost(2)*Mpi*Xc;
    Fc = (xx2.*xc)*One'- (xx2*xc');
    Fc = (Fc+Fc')/2 ;
    F=F-Fc.*P;
    Mis(1,2)=0;
    for i=1:dim
        if xx2(i)>0
            Mis(1,2)=Mis(1,2)+Cost(2)*((xx2(i)-meanX2P).^2)/2;
        else
            Mis(1,2)=Mis(1,2)+Cost(2)*((xx2(i)-meanX2M).^2)/2;
        end
    end
else
    Mis(1,2)=Cost(2);
end

if Dcost(3)>0
    
    [Xm,meanXmP]=constr_matrix_m(x_fix,x2');
    xm = Cost(3)*Mpi*Xm;
    
    Fm = (x2.*xm)*One'- (x2*xm');
    Fm = (Fm+Fm')/2 ;
    F=F+Fm.*P;
    
    Mis(1,3)=Cost(3)*meanXmP^2;
    
else
    Mis(1,3)=Cost(3);
end

% projection;
F =  F - diag(diag(F));
G = -F;

ng = norm(G,'fro');

% Norm conservation
alpha = -real(trace(E0'*G));
G = G + alpha*E0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% release active constraints if antigradient component G(i,j) is positive %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
    for j=1:i-1
        if Z(i,j)==0
            if G(i,j)> 0
                Z(i,j)=1;
                Z(j,i)=1;
            end
        end
    end
end

% projection on active constraints
G = G.*Z;
% Normalization of the gradient
G = G/norm(G,'fro');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  R.H.S. of the ODE          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Edot = G;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda(1)=lam2;
mis(1)=sum(Mis(1,:));

if store==1
    sol(:,:,1) = E0;
end
grid(1)=t0;
lamold=lam2;
misold=mis(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while t0<tf
    h
    %   starting stepsize
    h = h0;
    if h < hmin
        %disp('minimal stepsize')
        return
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   E0 is the current perturbation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    Eh = E0 + h*Edot;
    Eh=Eh/norm(Eh,'fro');
    Eh=check_B02(W,E0,h,Edot,threshold,epsil);

    %%%%%%%%%%%%%%%%%%%%%%%
E=Eh;norm(E,'fro')

B=W+epsil*E



    if (B >= -facth*threshold) % some elements may be fixed to -threshold
        % 
    else
        min(min(B))
        norm(E,'fro')
        disp('W+epsil*E is not positive after correction.. restart')
        return
    end
    

    %------------------------
    L = Lap(B);
    [lamsort, X] = compeigvecs(L, Wtype, 0, zeros(n), 3);
    lam2=lamsort(2);x2=X(:,2);
    %     if X(:,2)>0
    % x2=X(:,1);lam2=lamsort(1);
    % elseif X(:,2)>0
    % x2=X(:,1);lam2=lamsort(1);
    % else
    % x2=X(:,2);lam2=lamsort(2);
    % end
   

    if Dcost(1)>0
        Mist(1,1)=Cost(1)*lam2;
    else
        Mist(1,1)=Cost(1);
    end
    if Dcost(2)>0
        XX2=sum((X(:,2)-mean(X(:,2))).^2);
        XX1=sum((X(:,1)-mean(X(:,1))).^2);
        if XX1>XX2
            xx2=X(:,1);
        else
            xx2=X(:,2);
        end
        [Xc,meanX2P,meanX2M]=constr_matrix2(xx2);
        
        Mist(1,2)=0;
        for i=1:dim
            if xx2(i)>0
                Mist(1,2)=Mist(1,2)+Cost(2)*((xx2(i)-meanX2P).^2)/2;
            else
                Mist(1,2)=Mist(1,2)+Cost(2)*((xx2(i)-meanX2M).^2)/2;
            end
        end
    else
        Mist(1,2)=Cost(2);
    end
    if Dcost(3)>0
        
        [Xm,meanXmP]=constr_matrix_m(x_fix,x2');
        
        Mist(1,3)=Cost(3)*meanXmP^2;
        
    else
        Mist(1,3)=Cost(3);
    end
    
    mistemp=sum(Mist);
    misold;
    
    
    %   Check monotonicity
    if ((mistemp)<(misold))
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %          Step is accepted
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %          Store solution
        if store==1
            sol(:,:,ns+1)=E;
        end
        grid(ns+1)=t0+h;
        lambda(ns+1)=lam2;
        Mis(ns+1,:)=Mist;
        mis(ns+1)=sum(Mis(ns+1,:));
        ns=ns+1;
        lamold=lam2;
        misold=mis(ns);
        t0=t0+h;
        E0=E;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %          New vectorfield                     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [lamsort, X] = compeigvecs(L, Wtype, 0, zeros(n), 3);
        x2=X(:,2);lam2=lamsort(2);
        %     if X(:,2)>0
        % x2=X(:,1);lam2=lamsort(1);
        % elseif X(:,2)>0
        % x2=X(:,1);lam2=lamsort(1);
        % else
        % x2=X(:,2);lam2=lamsort(2);
        % end
        Mpi=pseudoinv(L,X,lam2,2);%pinv(L-lam2*eye(n));%Mpi=(Mpi+Mpi')/2;
        F=zeros(n);
        
        if Dcost(1)>0
            
            F2 = x2.^2*One';
            F2 = (F2+F2')/2 - (x2*x2');
            F = Cost(1)*(F2).*P;
            Mist(1,1)=Cost(1)*lam2;
        else
            Mist(1,1)=Cost(1);
        end
        if Dcost(2)>0
            if X(:,2)>0
                xx2=X(:,1);
            elseif X(:,2)<0
                xx2=X(:,1);
            else
                xx2=X(:,2);
            end
            [Xc,meanX2P,meanX2M]=constr_matrix2(xx2);
            xc = Cost(2)*Mpi*Xc;%(Mpi,x2I,x2);
            
            Fc = (xx2.*xc)*One'- (xx2*xc');
            Fc = (Fc+Fc')/2 ;
            F=F-Fc.*P;
            
            Mist(1,2)=0;
            for i=1:dim
                if xx2(i)>0
                    Mist(1,2)=Mist(1,2)+Cost(2)*((xx2(i)-meanX2P).^2)/2;
                else
                    Mist(1,2)=Mist(1,2)+Cost(2)*((xx2(i)-meanX2M).^2)/2;
                end
            end
            %Mist(1,2)=sum((xc-x2).^2)/2;
        else
            Mist(1,2)=Cost(2);
        end
        if Dcost(3)>0
            
            [Xm,meanXmP]=constr_matrix_m(x_fix,x2');
            xm = Cost(3)*Mpi*Xm;%(Mpi,x2I,x2);
            
            Fm = (x2.*xm)*One'- (x2*xm');
            Fm = (Fm+Fm')/2 ;
            F=F+Fm.*P;
            
            Mist(1,3)=Cost(3)*meanXmP^2;
            
        else
            Mist(1,3)=Cost(3);
        end
        
        
        F =  F - diag(diag(F));
        G = -F;
        % Norm conservation
        alpha = -real(trace(E0'*G));
        G = G + alpha*E0;
        % release active constraints if gradient component G(i,j) is positive
        for i=1:n
            for j=1:i-1
                if Z(i,j)==0
                    if G(i,j)> Edotmin % active constraint can be disactivated
                        Z(i,j)=1;
                        Z(j,i)=1;
                    end
                end
            end
        end
        %          Projection on active constraints
        G = G.*Z;
        %          Normalization of the gradient
        
        G = G/norm(G,'fro');
        
        
        %          New R.H.S.
        Edot = G;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if reject==0
            h0=h*gamma;
        else
            h0=h;
        end
        reject=0;
    else
        %           disp('rejected step')
        %           pause
        h0=h/gamma;
        reject=1;
    end
    
if misold-mistemp<=10*10^-5*h*misold+(10^-5)/100
    disp('fine')
return
end

    if G<threshold
        return
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Main loop                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if X(:,2)>0
%  xx2=X(:,1);
%  elseif X(:,2)<0
%  xx2=X(:,1);
%  else
%  xx2=X(:,2);
%  end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Laplacian of the graph           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L=Lap(M)

n=length(M);
One=ones(n,1);
L=diag(M*One)-M;

return
end
