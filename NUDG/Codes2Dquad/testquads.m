Globals2D;
N=8;
%[Nv,VX, VY, K, EToV,BCType,node,edge] = readmsh('mysquare.msh');
[Nv,VX, VY, K, EToV,BCType,node,edge] = readmsh('mysquare_refine.msh');
StartUp2Dquad;
figure(1);
plotmesh2dquad;
%hold on;
%plot(x,y,'.');

%test an integral
f= sin(pi*x/2).*sin(pi*y/2);
myint = dgintquad(f,V,J);
err = myint-4/pi^2


%BCType(BCType==Wall) = Dirichlet;
BuildBCMaps2D;
%try homogenous poisson problem
Lap = -PoissonIPDG2Dquad;

OP = [Lap ones(Np*K,1);
      ones(1,Np*K) 0];

uexact = cos(pi*x).*cos(pi*y);
f = -2*pi^2*cos(pi*x).*cos(pi*y);
%uexact = sin(pi*x).*sin(pi*y);
%f = -2*pi^2*sin(pi*x).*sin(pi*y);
rhs = MassMatrix*(J.*f);

u=OP\[rhs(:);0];
u = reshape(u(1:end-1),Np,K);

%u = Lap\rhs(:);
%u = reshape(u,Np,K);

figure(2);
subplot(2,1,1);
pf2dquad(N,x,y,u); shading interp; colorbar;
%subplt
% set up Dirichlet boundary conditions
% uD = zeros(Nfp*Nfaces, K);
% uD(mapD) = sin(pi*Fx(mapD)).*sin(pi*Fy(mapD));
% 
% % set up Neumann boundary conditions
% qN = zeros(Nfp*Nfaces, K);
% qN(mapN) = nx(mapN).*(pi*cos(pi*Fx(mapN)).*sin(pi*Fy(mapN))) + ...
%            ny(mapN).*(pi*sin(pi*Fx(mapN)).*cos(pi*Fy(mapN))) ;
% 
% % evaluate boundary condition contribution to rhs
% Aqbc = PoissonIPDGbc2D(uD, qN);

%rhs = MassMatrix*(J.*f) + Aqbc;

% for j=1:K
%     amp=0.005;
%     text(0.5*(VX(EToV(j,1))+VX(EToV(j,2)) )+randn(1,1)*amp,0.5*(VY(EToV(j,1))+VY(EToV(j,2)))+randn(1,1)*amp,'1');
%     text(0.5*(VX(EToV(j,2))+VX(EToV(j,3)) )+randn(1,1)*amp,0.5*(VY(EToV(j,2))+VY(EToV(j,3)))+randn(1,1)*amp,'2');
%     text(0.5*(VX(EToV(j,3))+VX(EToV(j,4)) )+randn(1,1)*amp,0.5*(VY(EToV(j,3))+VY(EToV(j,4)))+randn(1,1)*amp,'3');
%     text(0.5*(VX(EToV(j,4))+VX(EToV(j,1)) )+randn(1,1)*amp,0.5*(VY(EToV(j,4))+VY(EToV(j,1)))+randn(1,1)*amp,'4');
%     %pause;
% end
%hold off;

% figure(2);
% for k=1:K
%      plot(VX(EToV(k,[1 2 3 4 1])),VY(EToV(k,[1 2 3 4 1])),'.-k');
%      hold on;
%      plot(VX(EToV(EToE(k,1),[1 2 3 4 1])),VY(EToV(EToE(k,1),[1 2 3 4 1])),'.-r');
%      plot(VX(EToV(EToE(k,2),[1 2 3 4 1])),VY(EToV(EToE(k,2),[1 2 3 4 1])),'.-b');
%      plot(VX(EToV(EToE(k,3),[1 2 3 4 1])),VY(EToV(EToE(k,3),[1 2 3 4 1])),'.-g');
%      plot(VX(EToV(EToE(k,4),[1 2 3 4 1])),VY(EToV(EToE(k,4),[1 2 3 4 1])),'.-c');
%      axis([-.1 1.1 -.1 1.1]);
%      hold off;
%      drawnow;
%      
%      pause;
% end


