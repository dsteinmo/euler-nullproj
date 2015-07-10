function [OPDx,OPDy,MM] = DxDySEM2D(semtol)

Globals2D;

% build local volume mass matrix
MassMatrix = invV'*invV;

% build DG derivative matrices
MM  = zeros(K*Np*Np, 3);  OP = zeros(K*Np*Np*(1+Nfaces), 3);  

[semNodes,i,j] = uniquenodes(x,y,'rows');

DOF = length(i);

DG2SEM_map = reshape(j, Np, K);

% keyboard;
% global node numbering
entries = (1:Np*Np)'; entriesMM = (1:Np*Np)'; 
for k1=1:K 
  rows1 = DG2SEM_map(:,k1)*ones(1,Np);
  cols1 = rows1';

  % Build local operators  
  Dx = rx(1,k1)*Dr + sx(1,k1)*Ds;   Dy = ry(1,k1)*Dr + sy(1,k1)*Ds;

  OP11 = J(1,k1)*(MassMatrix*Dx);  %volume contribution to the laplacian.
                                                           %...transpose?
                                                           %Yes, because operator is stored in
  
  OP22 = J(1,k1)*(MassMatrix*Dy);
              
  %column vector temporarily.
  %keyboard;
  OPDx(entries(:), :)   = [rows1(:), cols1(:), OP11(:)];
  OPDy(entries(:), :)   = [rows1(:), cols1(:), OP22(:)];
  MM(entriesMM(:), :) = [rows1(:), cols1(:), J(1,k1)*MassMatrix(:)];
  entries = entries + Np*Np; entriesMM = entriesMM + Np*Np;
end  

OPDx   =   OPDx(1:max(entries)  -Np*Np,:);  OPDx   = myspconvert(OPDx, DOF, DOF, 1e-15);
OPDy   =   OPDy(1:max(entries)  -Np*Np,:);  OPDy   = myspconvert(OPDy, DOF, DOF, 1e-15);
MM   =   MM(1:max(entriesMM)-Np*Np,:);  MM   = myspconvert(MM, DOF, DOF, 1e-15);

return
