%BUG: Something is causing Fscale to be negative when it shouldn't
%e.g., if y is negative then it seems that Fscale is negative too. Find out why.

% Purpose : Setup script, building operators, grid, metric, and connectivity tables.
% Definition of constants
Nfp = N+1; Np = (N+1)*(N+1); Nfaces=4; 
NODETOL = 1e-7;

% Compute nodal set
[r,s] = Nodes2Dquad(N);

% Build reference element matrices
V = Vandermonde2Dquad(N,r,s);

invV = inv(V);
MassMatrix = invV'*invV;

[Dr,Ds] = Dmatrices2Dquad(N, r, s, V);

%I'm this far.

% build coordinates of all the nodes
v1 = EToV(:,1)'; v2 = EToV(:,2)'; v3 = EToV(:,3)'; v4 = EToV(:,4)';

%note that this is a _rectangular_ mapping (not general quad)
%could use something more general that uses all 4 vertices.
  x = 0.5*(r+1)*(VX(v2)-VX(v1))+repmat(VX(v1),length(r),1);
  y = 0.5*(s+1)*(VY(v3)-VY(v2))+repmat(VY(v2),length(s),1);

%triangle stuff:
%va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
%x = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));
%y = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));

% find all the nodes that lie on each edge
fmask1   = find( abs(s+1) < NODETOL)'; 
fmask2   = find( abs(r-1) < NODETOL)';
fmask3   = find( abs(s-1) < NODETOL)';
fmask4   = find( abs(r+1) < NODETOL)';
Fmask  = [fmask1;fmask2;fmask3;fmask4]';
Fx = x(Fmask(:), :); Fy = y(Fmask(:), :);

% Create surface integral terms
%%TODO: this.
LIFT = Lift2Dquad();

% calculate geometric factors. (this should work)
[rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);

% calculate geometric factors
[nx, ny, sJ] = Normals2Dquad(); %passed 1 test.
Fscale = sJ./(J(Fmask,:));

% Build connectivity matrix (done. passes one test also).
[EToE, EToF] = tiConnect2Dquad(EToV);

% Build connectivity maps
BuildMaps2D;
%BuildMaps2Dquad;  %turns out didn't need a new version.

% Compute weak operators (could be done in preprocessing to save time)
[Vr, Vs] = GradVandermonde2Dquad(N, r, s);
Drw = (V*Vr')/(V*V'); Dsw = (V*Vs')/(V*V');
