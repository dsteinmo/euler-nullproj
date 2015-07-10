function [nx, ny, sJ] = Normals2Dquad()

% function [nx, ny, sJ] = Normals2D()
% Purpose : Compute outward pointing normals at elements faces and surface Jacobians

Globals2D;

xr = Dr*x; yr = Dr*y; xs = Ds*x; ys = Ds*y; J = xr.*ys-xs.*yr;

% interpolate geometric factors to face nodes
fxr = xr(Fmask, :); fxs = xs(Fmask, :); fyr = yr(Fmask, :); fys = ys(Fmask, :);

% build normals
nx = zeros(Nfaces*Nfp, K); ny = zeros(Nfaces*Nfp, K);
fid1 = (1:Nfp)'; fid2 = (Nfp+1:2*Nfp)'; fid3 = (2*Nfp+1:3*Nfp)';
%new (for quads):
fid4 = (3*Nfp+1:4*Nfp)';
%should be able to use a big array and build this in a vectorized way.



%face 1: (top of quad?)
nx(fid1, :) =  fyr(fid1,:); ny(fid1,:) = fxr(fid1,:);

%face2: (left face?);
nx(fid2, :) = -fys(fid2,:); ny(fid2,:) = fxs(fid2,:);

%face 3: (bottom face)
nx(fid3, :) = fyr(fid3,:); ny(fid3,:)  =-fxr(fid3,:);

%face 4: (right face)
nx(fid4, :) = fys(fid4,:); ny(fid4,:)  = fxs(fid4,:);

%stuff seems backwards, so...:
nx = -nx;
ny = -ny;

%triangle shite:
% face 1
%nx(fid1, :) =  fyr(fid1, :); ny(fid1, :) = -fxr(fid1, :);

% face 2
%nx(fid2, :) =  fys(fid2, :)-fyr(fid2, :); ny(fid2, :) = -fxs(fid2, :)+fxr(fid2, :);

% face 3
%nx(fid3, :) = -fys(fid3, :); ny(fid3, :) =  fxs(fid3, :);

% normalise
sJ = sqrt(nx.*nx+ny.*ny); nx = nx./sJ; ny = ny./sJ;
return;
