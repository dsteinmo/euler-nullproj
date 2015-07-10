%get a mesh that doesn't have the stupid corners

clear;
close all;

% Driver script for solving the 2D Euler equations
Globals2D;

% Order of polynomials used for approximation 
N = 3;

% Read in Mesh
%load('../goodbox.mat');
load('../newxzplane2.mat');

%[Nv, VX, VY, K, EToV, BCType] = MeshReaderGambitBC2D('squarereg.neu');
%re-scale to unit square
% VX=(VX+1)/2;
% VY=(VY+1)/2;
% 
% % set up nodes and basic operations
% StartUp2D;
% 
% %h-refine
% refineflag= ones(K,1);
% Hrefine2D(refineflag);

errtol = 1e-11;

BCSolution = @NHSWMBC2D;

g=9.81;

% set up nodes and basic operations
StartUp2D;

% turn Neuman into walls
ids = find(BCType==7); 
BCType(ids) = Wall;

BuildBCMaps2D;

projectionMatrix = buildProjectionMatrix(N,r,s,V);

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
k1=66;
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

x0 = repmat(0.5*(xb+xc),Np,1);
y0 = repmat(0.5*(yb+yc),Np,1);
%rotate to (r,s) coordinates
ur = Ainv11.*(u-x0) + Ainv12.*(v-y0); 
us = Ainv21.*(u-x0) + Ainv22.*(v-y0); 

figure(1);
subplot(2,1,1);
PlotField2D_1tri(N,x(:,k1),y(:,k1),u); colorbar;
subplot(2,1,2);
PlotField2D_1tri(N,x(:,k1),y(:,k1),v); colorbar;
drawnow;
figure(2);
subplot(2,1,1);
PlotField2D_1tri(N,r,s,ur); colorbar;
subplot(2,1,2);
PlotField2D_1tri(N,r,s,us); colorbar;
drawnow;


[urproj,usproj] = divFreeProjection(ur,us,projectionMatrix);
figure(3);
subplot(2,1,1);
PlotField2D_1tri(N,r,s,urproj); colorbar;
subplot(2,1,2);
PlotField2D_1tri(N,r,s,usproj); colorbar;
drawnow;

%rotate back to (x,y) coordinates
uproj = A11.*urproj + A12.*usproj + x0;
vproj = A21.*urproj + A22.*usproj + y0;

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

u = -x(:,k1);
v = y(:,k1);

ur = Ainv11.*(u-x0) + Ainv12.*(v-y0); 
us = Ainv21.*(u-x0) + Ainv22.*(v-y0); 

[urproj,usproj] = divFreeProjection(ur,us,projectionMatrix);

uproj = A11.*urproj + A12.*usproj + x0;
vproj = A21.*urproj + A22.*usproj + y0;

err = (norm(u-uproj,2) + norm(v-vproj,2))/Np
if err < errtol
    disp('test 4 PASSED');
else
    disp('test 4 FAILED');
end

%%OR: (more expensive storage-wise for global projections)
%%Build velocity basis tailor-made for this triangle:
%%might have to do this for curvilinear elements.
%%(have to store MassMat and upsi/vpsi for each curved triangle.)

[divFreeMassMat,upsi,vpsi] = buildDivFreeMassMatrix(N,r,s,V,J(:,k1),rx(:,k1),ry(:,k1),sx(:,k1),sy(:,k1));
[uproj,vproj] = divfreeproj(u,v,upsi,vpsi,divFreeMassMat,V,J(:,k1));

err = (norm(u-uproj,2) + norm(v-vproj,2))/Np

if err < errtol
    disp('test 5 PASSED');
else
    disp('test 5 FAILED');
end

%test 6: a global solid-body vortex

%build coordinate rotation information on all elements.
A11 = zeros(Np,K); A12 = zeros(Np,K);
A21 = zeros(Np,K); A22 = zeros(Np,K);
Ainv11 = zeros(Np,K); Ainv12 = zeros(Np,K);
Ainv21 = zeros(Np,K); Ainv22 = zeros(Np,K);
x0 = zeros(Np,K); y0 = zeros(Np,K);
%ten extra fields to store, altogether.

disp('building global coordinate rotation data...');
for k1=1:K
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
    
    A11(:,k1) = repmat(Amat(1,1),Np,1); A12(:,k1) = repmat(Amat(1,2),Np,1);
    A21(:,k1) = repmat(Amat(2,1),Np,1); A22(:,k1) = repmat(Amat(2,2),Np,1);
    Ainv11(:,k1) = repmat(Ainv(1,1),Np,1); Ainv12(:,k1) = repmat(Ainv(1,2),Np,1);
    Ainv21(:,k1) = repmat(Ainv(2,1),Np,1); Ainv22(:,k1) = repmat(Ainv(2,2),Np,1);

    x0(:,k1) = repmat(0.5*(xb+xc),Np,1);
    y0(:,k1) = repmat(0.5*(yb+yc),Np,1);
end
clear Amat Ainv xa xb xc ya yb yc
disp('done.');

rotData = buildCoordRotationData; %rotData is a struct.

u = -x;
v = y;

% %rotate (u,v) to (ur,us) coordinates
% ur = rotData.Ainv11.*(u-rotData.x0) + rotData.Ainv12.*(v-rotData.y0); 
% us = rotData.Ainv21.*(u-rotData.x0) + rotData.Ainv22.*(v-rotData.y0); 
% %project on standard element
% [urproj,usproj] = divFreeProjection(ur,us,projectionMatrix);
% %rotate back
% uproj = rotData.A11.*urproj + rotData.A12.*usproj + rotData.x0;
% vproj = rotData.A21.*urproj + rotData.A22.*usproj + rotData.y0;

[uproj,vproj] = globalDivFreeProjection(u,v,projectionMatrix,rotData);

err = (norm(u-uproj,2) + norm(v-vproj,2))/(2*Np*K)

if err < errtol
    disp('test 6 (global projection) PASSED');
else
    disp('test 6 (global projection) FAILED');
end

%test 7: a totally divergent field, divergence should be nullified
u = x.^2;
v = y.^2;

[uproj,vproj] = globalDivFreeProjection(u,v,projectionMatrix,rotData);

err = norm(Div2D(uproj,vproj),2)/(Np*K)

if err < errtol
    disp('test 7 (global projection) PASSED');
else
    disp('test 7 (global projection) FAILED');
end

%test 8: a combo  field, divergence should be nullified
u = -y + 1e-2*x.^2;
v =  x + 1e-2*y.^2;

[uproj,vproj] = globalDivFreeProjection(u,v,projectionMatrix,rotData);

err = norm(Div2D(uproj,vproj),2)/(Np*K)

if err < errtol
    disp('test 8 (global projection) PASSED');
else
    disp('test 8 (global projection) FAILED');
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