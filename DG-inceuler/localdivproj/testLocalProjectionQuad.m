%get a mesh that doesn't have the stupid corners

clear;
close all;

addpath ~/work/NUDG/Codes2Dquad/

% Driver script for solving the 2D Euler equations
Globals2D;

% Order of polynomials used for approximation 
N = 4;
[Nv,VX, VY, K, EToV,BCType,node,edge] = readmsh('mysquare.msh');

errtol = 1e-11;


% set up nodes and basic operations
StartUp2Dquad;

% turn Neuman into walls
ids = find(BCType==7); 
BCType(ids) = Wall;

BuildBCMaps2D;

projectionMatrix = buildProjectionMatrixquad(N,r,s,V);

%test 1 (solid-body vortex, projection should keep it intact):
u = -s;
v = r;

[uproj,vproj] = divFreeProjection(u,v,projectionMatrix);

err = norm(u-uproj,2) + norm(v-vproj,2);
if err < errtol
    disp('test 1 PASSED');
else
    disp('test 1 FAILED');
end

%test 2: a totally divergent field, divergence be nullified.
u=r.^2;
v=s.^2;

[uproj,vproj] = divFreeProjection(u,v,projectionMatrix);

err = Dr*uproj+Ds*vproj;
if err < errtol
    disp('test 2 PASSED');
else
    disp('test 2 FAILED');
end


%test3: global field
rotData = buildCoordRotationData_quad; %rotData is a struct.

u = -x;
v = y;


[uproj,vproj] = globalDivFreeProjection(u,v,projectionMatrix,rotData);

err = (norm(u-uproj,2) + norm(v-vproj,2))/(2*Np*K);

if err < errtol
    disp('test 3 (global projection) PASSED');
else
    disp('test 3 (global projection) FAILED');
end

%test 4: a totally divergent field, divergence should be nullified
u = x.^2;
v = y.^2;

[uproj,vproj] = globalDivFreeProjection(u,v,projectionMatrix,rotData);

err = norm(Div2D(uproj,vproj),2)/(Np*K);

if err < errtol
    disp('test 4 (global projection) PASSED');
else
    disp('test 4 (global projection) FAILED');
end

%test 5: a combo  field, divergence should be nullified
u = -y + 1e-2*x.^2;
v =  x + 1e-2*y.^2;

[uproj,vproj] = globalDivFreeProjection(u,v,projectionMatrix,rotData);

err = norm(Div2D(uproj,vproj),2)/(Np*K);

if err < errtol
    disp('test 5 (global projection) PASSED');
else
    disp('test 5 (global projection) FAILED');
end


[efuns,d] = eig(projectionMatrix);
d = real(diag(d));
[d,ind] = sort(d,'descend');
efuns = efuns(:,ind);

% for jj=1:length(efuns)
%     uefun = efuns(1:end/2,jj);
%     vefun = efuns(end/2+1:end,jj);
%     figure(1); clf;
%     subplot(2,1,1);
%     PlotField2D_1tri(N,r,s,uefun); colorbar;
%     title(['\lambda=' num2str(d(jj)) '. inf norm of div u = ' num2str(norm(Dr*uefun + Ds*vefun,inf))]);
%     subplot(2,1,2);
%     PlotField2D_1tri(N,r,s,vefun); colorbar;
%     drawnow;
%     pause;
% end