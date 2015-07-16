function [OP] = build_c1x_misfit()


Globals2D;

% build local face matrices
massEdge = zeros(Np,Np,Nfaces);
Fm = Fmask(:,1); faceR = r(Fm); 
V1D = Vandermonde1D(N, faceR);  identEdge(Fm,Fm,1) = eye(size(V1D)); massEdge(Fm,Fm,1) = inv(V1D*V1D');
Fm = Fmask(:,2); faceS = s(Fm); 
V1D = Vandermonde1D(N, faceS);  identEdge(Fm,Fm,2) = eye(size(V1D)); massEdge(Fm,Fm,2) = inv(V1D*V1D');
Fm = Fmask(:,3); faceR = r(Fm); 
V1D = Vandermonde1D(N, faceR);  identEdge(Fm,Fm,3) = eye(size(V1D)); massEdge(Fm,Fm,3) = inv(V1D*V1D');
Fm = Fmask(:,4); faceS = s(Fm);
V1D = Vandermonde1D(N, faceS);  identEdge(Fm,Fm,4) = eye(size(V1D)); massEdge(Fm,Fm,4) = inv(V1D*V1D');

% build local volume mass matrix
MassMatrix = invV'*invV;

% build DG derivative matrices
%note: here, 3 is a sparse-formatting thing. not dependent on Nfaces.
OP = zeros(K*Np*Np*(1+Nfaces), 3);  

% global node numbering
entries = (1:Np*Np)'; 
for k1=1:K 
  if(~mod(k1,1000)) k1, end;
  rows1 = ((k1-1)*Np+1:k1*Np)'*ones(1,Np); cols1 = rows1';

  % Build local operators  
  Dx = rx(1,k1)*Dr + sx(1,k1)*Ds;   Dy = ry(1,k1)*Dr + sy(1,k1)*Ds;

  OP11 = zeros(size(Dx));
  OP12 = zeros(size(Dx));
  % Build element-to-element parts of operator
  for f1=1:Nfaces
    k2 = EToE(k1,f1); f2 = EToF(k1,f1); 

    rows2 = ((k2-1)*Np+1:k2*Np)'*ones(1,Np); cols2 = rows2';
    
    fidM  = (k1-1)*Nfp*Nfaces + (f1-1)*Nfp + (1:Nfp);
    vidM = vmapM(fidM); Fm1 = mod(vidM-1,Np)+1;
    vidP = vmapP(fidM); Fm2 = mod(vidP-1,Np)+1;
    
    id = 1+(f1-1)*Nfp + (k1-1)*Nfp*Nfaces;
    lnx = nx(id);  lny = ny(id); lsJ = sJ(id);
    hinv = max(Fscale(id), Fscale(1+(f2-1)*Nfp, k2));    

    Dx2 = rx(1,k2)*Dr + sx(1,k2)*Ds;   Dy2 = ry(1,k2)*Dr + sy(1,k2)*Ds;
    
    Dn1 = lnx*Dx  + lny*Dy ;
    Dn2 = lnx*Dx2 + lny*Dy2;

    mmE = lsJ*massEdge(:,:,f1);

    ident = identEdge(:,:,f1);

    switch(BCType(k1,f1))
      case {Dirichlet}
	     % nada
      case {Neuman}
        % nada 
      case {Wall}
        % nada   
      otherwise
	% interior face variational terms
	OP11        = OP11 + ( mmE*Dx - Dx'*mmE );

	OP12 = zeros(Np);
	OP12(Fm1,:) = OP12(Fm1,:) - (      mmE(Fm1,Fm1)*Dx2(Fm2,:) );  %corresponds to line 71, term 2
	OP12(:,Fm2) = OP12(:,Fm2) - (-Dx'*mmE(:, Fm1) );              %corresponds to line 71, term 3

	OP(entries(:), :) = [rows1(:), cols2(:), OP12(:)];

	entries = entries + Np*Np;

    end 
  end      
  OP(entries(:), :)   = [rows1(:), cols1(:), OP11(:)];
  entries = entries + Np*Np;

end  

OP   =   OP(1:max(entries)  -Np*Np,:);  OP   = myspconvert(OP, Np*K, Np*K, 1e-15);
return
