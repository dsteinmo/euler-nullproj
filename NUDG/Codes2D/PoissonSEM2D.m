function [OP,MM] = PoissonSEM2D(semtol)

Globals2D;

% build local volume mass matrix
MassMatrix = invV'*invV;

% build DG derivative matrices
MM  = zeros(K*Np*Np, 3);  OP = zeros(K*Np*Np*(1+Nfaces), 3);  

[semNodes,i,j] = uniquenodes(x,y);

DOF = length(i);

DG2SEM_map = reshape(j, Np, K);

% global node numbering
entries = (1:Np*Np)'; entriesMM = (1:Np*Np)'; 
for k1=1:K 
  rows1 = DG2SEM_map(:,k1)*ones(1,Np);
  cols1 = rows1';

  % Build local operators  
  Dx = rx(1,k1)*Dr + sx(1,k1)*Ds;   Dy = ry(1,k1)*Dr + sy(1,k1)*Ds;

  OP11 = J(1,k1)*(Dx'*MassMatrix*Dx + Dy'*MassMatrix*Dy);  %volume contribution to the laplacian.
                                                           %...transpose?
                                                           %Yes, because operator is stored in
                
  %column vector temporarily.
  %keyboard;
  OP(entries(:), :)   = [rows1(:), cols1(:), OP11(:)];
  MM(entriesMM(:), :) = [rows1(:), cols1(:), J(1,k1)*MassMatrix(:)];
  entries = entries + Np*Np; entriesMM = entriesMM + Np*Np;
  %entries
  %pause;
end  

OP   =   OP(1:max(entries) - Np*Np,:);  OP   = myspconvert(OP, DOF, DOF, 0);
MM   =   MM(1:max(entriesMM)- Np*Np,:);  MM   = myspconvert(MM, DOF, DOF, 0);

return
