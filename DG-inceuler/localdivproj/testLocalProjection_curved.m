%get a mesh that doesn't have the stupid corners

clear;
close all;

% Driver script for solving the 2D Euler equations
Globals2D;

% Order of polynomials used for approximation 
N = 8;

% Read in Mesh
%load('../goodbox.mat');
%load('../newxzplane2.mat');
load('mesh2dcircle.mat');

BCSolution = @NHSWMBC2D;
StartUp2D;

dtheta = pi/80;  %specify spacing at the bdry.
theta  = (-pi:dtheta:(pi-dtheta))';

rmax = 636.6;
%outer nodes
node   = [rmax.*cos(theta) rmax.*sin(theta)];
node = [node; node(1,:)];

%arclength parameterization of boundary.
numsplines =1; %1 for noisl, 2 for with isl
xs = cell(2,1);
ys = cell(2,1);
ts = cell(2,1);

[xs{1},ys{1},ts{1}] = ParametricSpline(node(:,1),node(:,2));
%[xs{2},ys{2},ts{2}] = ParametricSpline(node2(:,1),node2(:,2));
%% adjust curved elements 
% curved = [];
for j=1:numsplines
    %[curved_tmp,blendfact_x_tmp,blendfact_y_tmp] = MakeCurvedEdges_derek(Wall,xs{j},ys{j},ts{j});
    curved_tmp = MakeCurvedEdges_derek(Wall,xs{j},ys{j},ts{j});
    curved=[curved;curved_tmp];
    BuildBCMaps2D;
end
curved = unique(sort(curved,'ascend'));
straight = setdiff(1:K,curved);


errtol = 1e-9;



g=9.81;

% set up nodes and basic operations
%StartUp2D;

% turn Neuman into walls
ids = find(BCType==7); 
BCType(ids) = Wall;

BuildBCMaps2D;

projectionMatrix = buildProjectionMatrix(N,r,s,V); %get reference triangle projection matrix.

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

%test 3: same deal, but different coordinates
k1=straight(5);
u=x(:,k1).^2;
v=y(:,k1).^2;

%rotate (u,v) to (ur,us) coordinates
va = EToV(k1,1)'; vb = EToV(k1,2)'; vc = EToV(k1,3)';
xa = VX(va); xb = VX(vb); xc = VX(vc);
ya = VY(va); yb = VY(vb); yc = VY(vc);

%x = 0.5*(-(r+s)*xa+(1+r)*xb+(1+s)*xc);
%y = 0.5*(-(r+s)*ya+(1+r)*yb+(1+s)*yc);
%or, 
%(x,y) = A *(r,s) + 0.5*(xb+xc,yb+yc)
%where:
Amat = 0.5*[(xb-xa) (xc-xa);
            (yb-ya) (yc-ya)];

Ainv =inv(Amat);
A11 = repmat(Amat(1,1),Np,1); A12 = repmat(Amat(1,2),Np,1);
A21 = repmat(Amat(2,1),Np,1); A22 = repmat(Amat(2,2),Np,1);
Ainv11 = repmat(Ainv(1,1),Np,1); Ainv12 = repmat(Ainv(1,2),Np,1);
Ainv21 = repmat(Ainv(2,1),Np,1); Ainv22 = repmat(Ainv(2,2),Np,1);

clear Ainv Amat

ur = Ainv11.*(u) + Ainv12.*(v); 
us = Ainv21.*(u) + Ainv22.*(v); 

% figure(1);
% subplot(2,1,1);
% PlotField2D_1tri(N,x(:,k1),y(:,k1),u); colorbar;
% subplot(2,1,2);
% PlotField2D_1tri(N,x(:,k1),y(:,k1),v); colorbar;
% drawnow;
% figure(2);
% subplot(2,1,1);
% PlotField2D_1tri(N,r,s,ur); colorbar;
% subplot(2,1,2);
% PlotField2D_1tri(N,r,s,us); colorbar;
% drawnow;


[urproj,usproj] = divFreeProjection(ur,us,projectionMatrix);
% figure(3);
% subplot(2,1,1);
% PlotField2D_1tri(N,r,s,urproj); colorbar;
% subplot(2,1,2);
% PlotField2D_1tri(N,r,s,usproj); colorbar;
% drawnow;

%rotate back to (x,y) coordinates
%uproj = A11.*urproj + A12.*usproj + x0;
%vproj = A21.*urproj + A22.*usproj + y0;

uproj = A11.*urproj + A12.*usproj;
vproj = A21.*urproj + A22.*usproj;

norm(Dr*urproj+Ds*usproj,2)

err = rx(:,k1).*(Dr*uproj) + sx(:,k1).*(Ds*uproj) ...
    + ry(:,k1).*(Dr*vproj) + sy(:,k1).*(Ds*vproj);

err = norm(err,2);


if err < errtol
    disp('test 3 PASSED');
else
    disp('test 3 FAILED');
end

%test 4 (solid-body vortex in different coords, projection should keep it intact):

u = -y(:,k1);
v = x(:,k1);

ur = Ainv11.*(u) + Ainv12.*(v); 
us = Ainv21.*(u) + Ainv22.*(v); 

[urproj,usproj] = divFreeProjection(ur,us,projectionMatrix);

uproj = A11.*urproj + A12.*usproj;
vproj = A21.*urproj + A22.*usproj;

err = (norm(u-uproj,2) + norm(v-vproj,2))/Np
if err < errtol
    disp('test 4 PASSED');
else
    disp('test 4 FAILED');
end


%%OR: (more expensive storage-wise for global projections)
%%Build velocity basis tailor-made for this triangle:

[divFreeMassMat,upsi,vpsi] = buildDivFreeMassMatrix(N,r,s,V,J(:,k1),rx(:,k1),ry(:,k1),sx(:,k1),sy(:,k1));
[uproj,vproj] = divfreeproj(u,v,upsi,vpsi,divFreeMassMat,V,J(:,k1));

err = (norm(u-uproj,2) + norm(v-vproj,2))/Np

if err < errtol
    disp('test 5 PASSED');
else
    disp('test 5 FAILED');
end

%for j=1:length(upsi)-1
%    disp(rx(:,k1).*(Dr*upsi(:,j))+sx(:,k1).*(Ds*upsi(:,j)) + ry(:,k1).*(Dr*vpsi(:,j))+sy(:,k1).*(Ds*vpsi(:,j)));
%end


disp(' ');
disp('test 6 (solid body vortex on curved element):');
%test 6:
%now try it for an actual curvilinear element:
%k2 = curved(1);
%k2 = k1;
k2 = curved(2);

%try with cubature rules
CubatureOrder = floor(2*(N+1)*3/2);  %orig
cub = CubatureVolumeMesh2D(CubatureOrder);

[cinfo] = BuildCurvedOPS2D(CubatureOrder); %this gives diff. matrices.

%curved(2)

%u = x(:,k2).^2/max(x(:))^2;
%v = y(:,k2).^2/max(y(:))^2;

u = -y(:,k2);
v = x(:,k2);

ux = cinfo(2).Dx*u;
vy = cinfo(2).Dy*v;


[divFreeMassMat,upsi,vpsi] = buildDivFreeMassMatrix_curved(cinfo(2).Dx,cinfo(2).Dy,V,cub.V,cub.W(:,curved(2)));
[uproj,vproj] = divfreeproj_curved(u,v,upsi,vpsi,divFreeMassMat,cub.V,cub.W(:,curved(2)));

% [divFreeMassMat,upsi,vpsi] = buildDivFreeMassMatrix_curved_monom(N,x(:,k2),y(:,k2),V,cub.V,cub.W(:,curved(2)));
% [uproj,vproj] = divfreeproj_curved(u,v,upsi,vpsi,divFreeMassMat,cub.V,cub.W(:,curved(2)));
% norm(cinfo(2).Dx*uproj + cinfo(2).Dy*vproj,2)

projectionMatrix_curved = buildProjectionMatrix_curved(cinfo(2).Dx,cinfo(2).Dy,V,cub.V,cub.W(:,curved(2)));
[uproj2,vproj2] = divFreeProjection_curved(u,v,projectionMatrix_curved);

disp('difference between direct projection and using a projection matrix:');
err = (norm(uproj-uproj2,2) + norm(vproj-vproj2,2))/Np

disp('error in projected result (via proj matrix) vs. expected result:');
err = (norm(u-uproj2,2) + norm(v-vproj2,2))/Np

figure(1);
subplot(2,1,1);
PlotField2D_1tri(N,x(:,k2),y(:,k2),u-uproj);
subplot(2,1,2);
PlotField2D_1tri(N,x(:,k2),y(:,k2),v-vproj);
drawnow;


disp('test 7 (fully compressible field on curved element):');
u = x(:,k2);
v = y(:,k2);
[uproj,vproj] = divFreeProjection_curved(u,v,projectionMatrix_curved);

%divergence should be nullified
divuproj = cinfo(2).Dx*uproj + cinfo(2).Dy*vproj;
err = norm(divuproj,2)



disp('test 8 (global solid body vortex, curved & straight elements');

u = -x;
v = y;

%pre-processing:
rotDataStraight = buildCoordRotationDataStraight; %rotData is a struct.
projectionMatricesCurved = zeros(2*Np,2*Np,length(curved));
disp('building & storing curved projection matrices...');
for n=1:length(curved)
        k1 = curved(n);
        projectionMatricesCurved(:,:,n)= buildProjectionMatrix_curved(cinfo(n).Dx,cinfo(n).Dy,V,cub.V,cub.W(:,k1));
end
disp('done.');
projMatrixStraight = projectionMatrix;

uproj = zeros(Np,K); vproj= zeros(Np,K);


%%do projection on straight elements:
% [uproj(:,straight),vproj(:,straight)] = globalDivFreeProjection(u(:,straight),v(:,straight),projectionMatrix,rotDataStraight);
% %do projection on each curved element:
% for n=1:length(curved)
%     k1=curved(n);
%     [uproj(:,k1),vproj(:,k1)] = divFreeProjection_curved(u(:,k1),v(:,k1),projectionMatricesCurved(:,:,n));
% end
[uproj,vproj] = globalDivFreeProjection_curved(u,v,projMatrixStraight,rotDataStraight,projectionMatricesCurved);

% figure(1); clf;
% subplot(2,1,1);
% pf2d(N,x,y,u-uproj); shading flat; colorbar;
% subplot(2,1,2);
% pf2d(N,x,y,v-vproj); shading flat; colorbar;


disp('error in projected result vs. expected result:');
err = (norm(u-uproj,2) + norm(v-vproj,2))/Np



%test 8: a totally divergent field, divergence should be nullified
u = x.^2/1e6;
v = y.^2/1e6;

[uproj,vproj] = globalDivFreeProjection_curved(u,v,projMatrixStraight,rotDataStraight,projectionMatricesCurved);

err = norm(Div2D_curved(uproj,vproj),2)/(Np*K)

if err < errtol
    disp('test 8 (global projection) PASSED');
else
    disp('test 8 (global projection) FAILED');
end

figure(1); clf;
subplot(3,2,1);
pf2d(N,x,y,u); colorbar;
subplot(3,2,2);
pf2d(N,x,y,uproj); colorbar;
subplot(3,2,3);
pf2d(N,x,y,v); colorbar;
subplot(3,2,4);
pf2d(N,x,y,vproj); colorbar;
subplot(3,2,5);
pf2d(N,x,y,Div2D_curved(u,v)); colorbar;
subplot(3,2,6);
pf2d(N,x,y,Div2D_curved(uproj,vproj)); colorbar;



%test 9: a combo  field, divergence should be nullified
u = (-y + 1e-2*x.^2)/1e3;
v =  (x + 1e-2*y.^2)/1e3;

[uproj,vproj] = globalDivFreeProjection_curved(u,v,projMatrixStraight,rotDataStraight,projectionMatricesCurved);

err = norm(Div2D(uproj,vproj),2)/(Np*K)

if err < errtol
    disp('test 9 (global projection) PASSED');
else
    disp('test 9 (global projection) FAILED');
end

% [efuns,d] = eig(projectionMatrix);
% d = real(diag(d));
% [d,ind] = sort(d,'descend');
% efuns = efuns(:,ind);

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