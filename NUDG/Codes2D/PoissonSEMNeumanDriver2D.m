% Driver script for solving the 2D Poisson equation
% Note that if you use Neumann conditions the solution will have some
% unknown constant added to it (non-uniqueness of potential function)
% Spatially dependent diffusivity working. 

Globals2D;

% Polynomial order used for approximation 
N = 4;

% Read in Mesh
%[Nv, VX, VY, K, EToV] = MeshReaderGambitBC2D('circA01.neu'); %Dirichlet circle (no curvilinear elements)
load('kinneret_donut_mesh.mat');

% Initialize solver and construct grid and metric
StartUp2D;
BuildBCMaps2D;

for zzz=1:0
   Hrefine2D(ones(K,1))
   StartUp2D;
   BuildBCMaps2D;
end

% set up right hand side for homogeneous Poisson 
[A,M] = PoissonSEM2D();

% set up Neumann boundary conditions (this only kicks in if mesh supports it)
qN = zeros(Nfp*Nfaces, K);
%qN(mapW) = nx(mapW).*(pi*cos(pi*Fx(mapW)/2000).*sin(pi*Fy(mapW)/2000)) + ...
%           ny(mapW).*(pi*sin(pi*Fx(mapW)/2000).*cos(pi*Fy(mapW)/2000)) ;

%qN(mapW) = qN(mapW)/1000;

rmax = max(x(:))
qN(mapW) = ny(mapW).*Fy(mapW)/(rmax*10);

% evaluate boundary condition contribution to rhs
%Aqbc = PoissonSEMNeumanbc2D(qN);


rmax = max(x(:));
% set up right hand side forcing
force = (1/rmax)*cos(pi*x/2000).*cos(pi*y/2000);  %use this for a = 1


nodes = [x(:) y(:)];
[semNodes,i,j] = unique(nodes, 'rows');

DOF = length(i);

DG2SEM_map = reshape(j, Np, K);

u = zeros(DOF,1);

forcesem = force(i);
rhssem = M*forcesem;

rhssem = rhssem + Aqbc(i);

%'Border' matrix 'A' such that mean of u is zero

OP = [A ones(DOF,1);
      ones(1,DOF) 0;];

u = OP\[rhssem;0];

sigma = u(end);
usem = u(1:end-1);

disp(['sigma: ' num2str(sigma) ] )
disp(['mean u: ' num2str(mean(u(:))) ] )

%scatter back to DG grid.
u = u(j);
u = reshape(u, Np, K);

uexact = sin(pi*x).*sin(pi*y);  %exact solution
%uexact = cos(pi*x).*cos(pi*y);


 figure(1);
 pf2d(N,x,y,u); 
 colormap(darkjet);
 colorbar;
 title('SEM Solution');
 

[OPDx,OPDy,M] = DxDySEM2D();

lumpM = spdiags(sum(M,2), [0], length(M), length(M));

disp('taking dewivatives...');

f = exp(-(x/2000).^2 -(y/2000).^2);

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
subplot(2,1,1);
pf2d(N,x,y,ux); colorbar; colormap(darkjet);
subplot(2,1,2);
pf2d(N,x,y,uy); colorbar; colormap(darkjet);
 
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
