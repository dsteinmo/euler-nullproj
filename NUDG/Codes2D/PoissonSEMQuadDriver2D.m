% Driver script for solving the 2D Poisson equation
% Note that if you use Neumann conditions the solution will have some
% unknown constant added to it (non-uniqueness of potential function)
% Spatially dependent diffusivity working. 
clear;
close all;
Globals2D;
addpath '../';
setNUDGpaths;
addpath '/home/dsteinmo/work/NUDG/Codes2Dquad/';
addpath 'SEMEuler2D/consolidator';
% Polynomial order used for approximation 
N = 4;

% Note: There are some 'bad elements' in lower left corner of 'testquad.msh.'
% Need to make a better mesh, or fix my quad mapping.
% Read in Mesh %testquad.msh, simple_quad.msh
[Nv, VX, VY, K, EToV, BCType,node,edge] = readmsh('testquad.msh');
%[Nv, VX, VY, K, EToV] = MeshReaderGambitBC2D('circA01.neu'); %Dirichlet circle (no curvilinear elements)
%load('kinneret_donut_mesh.mat');

%EToV(3,:) = [EToV(3,2:4) EToV(3,1)];
% Initialize solver and construct grid and metric
%StartUp2Dquad;
%BCType = CorrectBCTablequad(EToV, find(abs(VY-min(y(:)))<NODETOL | abs(VY-max(y(:)))<NODETOL | abs(VX-min(x(:)))<NODETOL | abs(VX-max(x(:)))<NODETOL),Wall,K);

figure(1);

%hack to make adjust local node orderings such that they will work
%with my non-general mapping.
found=0;
for k=1:K
   if (abs(VX(EToV(k,1)) - VX(EToV(k,2)))<1e-3)
      EToV(k,:) = [EToV(k,2:4) EToV(k,1)];
      found = found + 1;
   end
   EToV(k,:) = EToV(k,4:-1:1);
end

StartUp2Dquad;
BuildBCMaps2D;

BCType = CorrectBCTablequad(EToV, find(abs(VY-min(y(:)))<NODETOL | abs(VY-max(y(:)))<NODETOL | abs(VX-min(x(:)))<NODETOL | abs(VX-max(x(:)))<NODETOL),Wall,K);

StartUp2Dquad;
BuildBCMaps2D;

%figure(1);
%plot(Fx(mapW),Fy(mapW),'.');
%drawnow;

% set up right hand side for homogeneous Poisson 
semtol = 9.9e-2; %5e-2 too big, 3e-2 too small
[A,M] = PoissonSEM2D(semtol);
%A=-A;
% set up Neumann boundary conditions (this only kicks in if mesh supports it)
qN = zeros(Nfp*Nfaces, K);
%qN(mapW) = nx(mapW).*(pi*cos(pi*Fx(mapW)/2000).*sin(pi*Fy(mapW)/2000)) + ...
%           ny(mapW).*(pi*sin(pi*Fx(mapW)/2000).*cos(pi*Fy(mapW)/2000)) ;

%qN(mapW) = qN(mapW)/1000;

%rmax = max(x(:))
%qN(mapW) = ny(mapW).*Fy(mapW)/(rmax*10);

% evaluate boundary condition contribution to rhs
%Aqbc = PoissonSEMNeumanbc2D(qN);

% set up right hand side forcing
%force = (1/rmax)*cos(pi*x/2000).*cos(pi*y/2000);  %use this for a = 1
force = -2*pi^2*sin(pi*x).*sin(pi*y);

nodes = [x(:) y(:)];

[semNodes,i,j] = uniquenodes(x,y);

DOF = length(i);

DG2SEM_map = reshape(j, Np, K);

forcesem = force(i);
rhssem = M*forcesem;

OP = [A ones(DOF,1);
      ones(1,DOF) 0;];

u = OP\[rhssem;0];

sigma = u(end)
usem = u(1:end-1);

%disp(['sigma: ' num2str(sigma) ] )
%disp(['mean u: ' num2str(mean(u(:))) ] )


%scatter back to DG grid.
u = usem(j);
u = reshape(u, Np, K);

uexact = cos(pi*x).*cos(pi*y);


 figure(1);
 pf2dquad(N,x,y,u); 
 colormap(darkjet);
 colorbar;
 title('SEM Solution');
 shading interp;


[OPDx,OPDy,M] = DxDySEM2D(semtol);
%[~,i,j] = unique([x(:) y(:)], 'rows');

lumpM = spdiags(sum(M,2), [0], length(M), length(M));

disp('taking dewivatives...');

f = exp(-(5*(x-.5)).^2 -(5*(y-.5)).^2);
%f = u;

disp('factorizing mass matrix...');
[ll,uu,pp,qq] = lu(M);
disp('done.');

%ux = M\(OPDx*f(i));
disp('calculating Dx...');
ux = qq*(uu\(ll\(pp*(OPDx*f(i)))));
ux = ux(j);
ux = reshape(ux,Np,K);
disp('done.');
%uy = M\(OPDy*f(i));
disp('Calculating Dy...');
uy = qq*(uu\(ll\(pp*(OPDy*f(i)))));
uy = uy(j);
uy = reshape(uy,Np,K);
disp('done.');

figure(2);
pf2dquad(N,x,y,ux); colorbar; colormap(darkjet);
shading interp;

figure(3);
pf2dquad(N,x,y,uy); colorbar; colormap(darkjet);
shading interp;
return;

uxx = qq*(uu\(ll\(pp*(OPDx*ux(i)))));
uxx = reshape(uxx(j),Np,K);

uyy = qq*(uu\(ll\(pp*(OPDy*uy(i)))));
uyy = reshape(uyy(j),Np,K);

figure(3);
subplot(2,1,1);
pf2dquad(N,x,y,uxx); colorbar; colormap(darkjet);
subplot(2,1,2);
pf2dquad(N,x,y,uyy); colorbar; colormap(darkjet);
 
 return;
% figure(2);
% PlotMesh2D;
% figure(3);
% PlotField2D(N,x,y,force); view([0 90]); colorbar;
% title('forcing');
% figure(4);
% PlotField2D(N,x,y,a); view([0 90]); colorbar;
% title('diffusivity');
% figure(5);
% PlotField2D(N,x,y,uexact); view([0 90]); colorbar;
% title('Exact Solution');
% error = norm(uexact(:)-u(:),2)/norm(uexact(:),2)


%look at eigenvalues
% numeigs=21;
% 
% [V,d]=eigs(A,numeigs,'SM');   
% %if using 'eig'
% % [V,d] = eig(full(A));
% %    d=diag(d);
% %    [d,inds] = sort(d,'descend');
% %    V = V(:,inds);
% %end if using 'eig
% phi = cell(numeigs,1);
% for jj=1:numeigs
%     phi{jj} = V(:,jj);
%     phi{jj} = reshape(phi{jj},Np,K);
%     figure(jj);
%     PlotField2D(N,x,y,phi{jj}); view([0 90]); colorbar;
% end
