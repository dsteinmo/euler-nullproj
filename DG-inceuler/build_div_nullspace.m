function [OP] = build_div_nullspace()


Globals2D;

% build local volume mass matrix
MassMatrix = invV'*invV;

% build DG derivative matrices
%note: here, 3 is a sparse-formatting thing. not dependent on Nfaces.
OP = zeros(K*Np*Np*(1+Nfaces), 3);  

% global node numbering
dof=Np*K

%loop over elements. 

count=0;
for k1=1:K 
  if(~mod(k1,1000)) k1, end;
  rows1 = ((k1-1)*Np+1:k1*Np)'*ones(1,Np); cols1 = rows1';

  % Build local operators  
  Dx = rx(1,k1)*Dr + sx(1,k1)*Ds;   Dy = ry(1,k1)*Dr + sy(1,k1)*Ds;

  Div = [Dx zeros(Np,Np);
       zeros(Np,Np) Dy];

  nullDiv = null(full(Div));
  %split up component-wise
  nullDiv_u = nullDiv(1:end/2,:);
  nullDiv_v = nullDiv(end/2+1:end,:);

  dims = size(nullDiv_u);
  width = dims(2);

  rows1 = rows1(:,1:width);
  cols1 = cols1(:,1:width);
   
  entries = count + (1:length(rows1(:)));
  count = count + (1:length(rows1(:)));
  OP(entries(:), :)   = [rows1(:), cols1(:), nullDiv_u(:)];
  
  entries = count + (1:length(rows1(:)));
  count = count + (1:length(rows1(:)));

  OP(entries(:), :) = [rows1(:)+dof, cols1(:)+dof, nullDiv_v(:)];
end  

OP   =   OP(1:max(entries) ,:);  OP   = myspconvert(OP, 2*Np*K, 2*Np*K, 1e-15);
return
