function [BigDx BigDy] = DxDyquad()

Globals2D;

% build local volume mass matrix
MassMatrix = invV'*invV;

% build DG derivative matrices
%note: here, 3 is a sparse-formatting thing. not dependent on Nfaces.
BigDx  = zeros(K*Np*Np, 3);  BigDy = zeros(K*Np*Np, 3);  

% global node numbering
entries = (1:Np*Np)'; entriesMM = (1:Np*Np)'; 
for k1=1:K 
  if(~mod(k1,1000)) k1, end;
  rows1 = ((k1-1)*Np+1:k1*Np)'*ones(1,Np); cols1 = rows1';

  % Build local operators  
  Dx = rx(1,k1)*Dr + sx(1,k1)*Ds;   Dy = ry(1,k1)*Dr + sy(1,k1)*Ds;

  BigDx(entries(:), :)   = [rows1(:), cols1(:), Dx(:)];
  BigDy(entries(:), :)   = [rows1(:), cols1(:), Dy(:)];
  entries = entries + Np*Np;
end  

BigDx   =   BigDx(1:max(entries)  -Np*Np,:); 

BigDy   =   BigDy(1:max(entries)  -Np*Np,:); 
BigDx   = myspconvert(BigDx, Np*K, Np*K, 1e-15);
BigDy   = myspconvert(BigDy, Np*K, Np*K, 1e-15);
return
