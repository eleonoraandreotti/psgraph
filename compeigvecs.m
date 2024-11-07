function [lamsort, X] = compeigvecs(A, Atype, epsil, E, m)
% get eigenvalues and eigenvectors of a real symmetric matrix
% get m smallest eigenvalues and correspoding eigenvectors X
% of B=A+epsil*E, with ||x||=1.
% if Atype == 1, A is dense so use dense matrix operations
% if Atype == 2, A is sparse matrix
% if Atype == 3, A is function handle



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Authors: Nicola Guglielmi and Christian Lubich
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    m = 1;
end

if Atype == 1 % A is dense, so compute all the eigenvalues
  if E==0
    [X, Lambda] = eig(A);
  else
    B = A + epsil*E;
    [X, Lambda] = eig(B);  
  end
else
  opts.disp = 0; % for eigs
  opts.tol = 1e-10;
  opts.isreal = 1; % for eigs, when Atype== 3
  opts.issym = 1;  % for eigs, when Atype== 3
% opts.v0 = X(:,1);  % possible input 
  eigscode = 'SA'; % want eigenvalue(s) with smallest magnitudes  
  if E==0
    [X, Lambda] = eigs(A, m, eigscode, opts); 
  else
    B = A + epsil*E;    
    %
    % compute eigenvectors
    %   Afun = @(p)aplusYYt(A, Atype, m, epsil, Y, alpha, p, PIJ);
    [X, Lambda] = eigs(B, m, eigscode, opts);
  end
end
lambda = diag(Lambda);
% sort eigenvalues: they are not sorted by eig (eigs), and although they
% are resorted
[lamsort, indx] = sort(lambda, 'ascend');
lamsort=lamsort(1:m);
X=X(:,indx(1:m));
