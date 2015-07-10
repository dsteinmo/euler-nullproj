%try a helmholtz decomposition


clear;
close all;

% Driver script for solving the 2D Euler equations
Globals2D;

% Order of polynomials used for approximation 
N = 4;

load('raobox.mat');

g=9.81;
H=max(y(:)); %'1'
L=max(x(:)); %'1'

BCSolution = @NHSWMBC2D;

StartUp2D;
% refineflag= ones(K,1);
% Hrefine2D(refineflag);
% StartUp2D;

% turn Neuman into walls
ids = find(BCType==7); 
BCType(ids) = Wall;

BuildBCMaps2D

%urot = -(y-.5).*exp(-10*((x-.5).^2+(y-.5).^2));
%vrot = (x-.5).*exp(-10*((x-.5).^2+(y-.5).^2));

urot = -(y-0.5).*exp(-50*((x-0.5).^2+(y-0.5).^2));
vrot = (x-0.5).*exp(-50*((x-0.5).^2+(y-0.5).^2));

udiv = (x-0.5).*exp(-50*(x-0.5).^2-(y-0.5).^2);
vdiv = (y-0.5).*exp(-50*(x-0.5).^2-(y-0.5).^2);

u = urot+1e-2*udiv;
v = vrot+1e-2*vdiv;

% figure(1);
% PlotField2D(N,x,y,u); view([0 90]); axis tight; shading interp; colorbar; drawnow;
% 
% figure(2);
% PlotField2D(N,x,y,v); view([0 90]); axis tight; shading interp; colorbar; drawnow;
% 
% %compute vorticity of our field (identically 2, same as 'rot' field)
% [dump1,dump2,vort] = Curl2D(u,v,[]);
% figure(3);
% PlotField2D(N,x,y,vort); view([0 90]); axis tight; shading interp; colorbar; drawnow;

%compute divergence of our field (identically 2, same as 'div' field)
div = Div2D(u,v);
figure(4);
PlotField2D(N,x,y,div); view([0 90]); axis tight; shading interp; colorbar; drawnow;

%build Neuman Laplacian
[OpNeu,MM] = PoissonIPDG2D_eig; %remember, this thing returns the positive definite operator, i.e., negative laplacian

numpot=10;
numstrm=75;

[efuns,d,convflag]=eigs(OpNeu,MM,numpot+1,'SM');   
d= diag(d);
oldd =d;
[lambda,inds] = sort(d,'ascend');
%[d,inds] = sort(d,'ascend');
efuns = efuns(:,inds);

%get rid of mode zero
lambda = lambda(2:end);
efuns = efuns(:,2:end);

phi = cell(numpot,1);  %Allocate cell array for homogenous Neuman eigenfns

for jj=1:numpot
    phitmp = efuns(:,jj);
    phi{jj} = reshape(phitmp,Np,K);
    %normalize eigenfunction
    alpha = dgint(phi{jj}.*phi{jj},V,J);
    phi{jj} = (sqrt(g*H^3*L^2)/sqrt(alpha)/sqrt(lambda(jj))).*phi{jj};
    
    %disp((lambda(jj)*dgint(phi{jj}.*phi{jj},V,J)-g*H^3*L^2));
end

% turn Neuman into walls
ids = BCType==Wall; 
BCType(ids) = Dirichlet;

BuildBCMaps2D;
%build dirichlet laplacian.
[OpDir,MM] = PoissonIPDG2D_eig;
ids = find(BCType==Dirichlet); 
BCType(ids) = Wall;

BuildBCMaps2D;

[efuns,d,convflag]=eigs(OpDir,MM,numstrm,'SM');   
if convflag ~= 0
    disp('eigenvalues didn''t converged! (normal laplacian operator)');
end

d= diag(d);
oldd =d;
[mu,inds] = sort(d,'ascend');
efuns = efuns(:,inds);


disp('checking normalization of psis...');
psi = cell(numstrm,1);  %Allocate cell array for homogenous Dirichlet eigenfns
for jj=1:numstrm
    psitmp = efuns(:,jj);
    psi{jj} = reshape(psitmp,Np,K); 
    
    %normalize eigenfunction
    alpha = dgint(psi{jj}.*psi{jj},V,J);
    psi{jj} = (sqrt(g*H^3*L^2)/sqrt(alpha)/sqrt(mu(jj))).*psi{jj};
    
    %disp(mu(jj)*dgint(psi{jj}.*psi{jj},V,J)-g*H^3*L^2);
end

%build solenoidal velocity basis
upsi = cell(numstrm,1);
vpsi = cell(numstrm,1);

dpsi = zeros(3*Nfp,K);
for jj=1:numstrm
    psitmp = psi{jj};
    [psix,psiy] = Grad2D(psitmp);
    dpsi(:) = psitmp(vmapM) - psitmp(vmapP);
    fluxpsix = nx.*(dpsi/2); fluxpsiy = ny.*(dpsi/2);
    
    %add surface integral contributions to volumetric gradient
    psix = psix - LIFT*(Fscale.*fluxpsix);
    psiy = psiy - LIFT*(Fscale.*fluxpsiy);
    
    upsi{jj} = -psiy;
    vpsi{jj} = psix;
    
end

% figure(1);
% PlotField2D(N,x,y,upsi{3}); view([0 90]); axis tight; shading interp; colorbar; drawnow;
% figure(2);
% PlotField2D(N,x,y,vpsi{3}); view([0 90]); axis tight; shading interp; colorbar; drawnow;

%project velocity field onto divergence-free basis functions
%i.e., compute Galerkin expansion coefficients and sum the
%truncated series in the same loop. (they're orthonormal)
c = zeros(numstrm,1);
uproj = zeros(Np,K); vproj = zeros(Np,K);
for jj=1:numstrm
    c(jj) = dgint(u.*upsi{jj} + v.*vpsi{jj},V,J)/(g*H^3*L^2);
    
    uproj = uproj + c(jj)*upsi{jj};
    vproj = vproj + c(jj)*vpsi{jj};
end

%form projected velocity

% for jj=1:numstrm
%     uproj = uproj + c(jj)*upsi{jj};
%     vproj = vproj + c(jj)*vpsi{jj};
% end

close all;

figure(1);
PlotField2D(N,x,y,uproj); view([0 90]); axis tight; shading interp; colorbar; drawnow;
figure(2);
PlotField2D(N,x,y,vproj); view([0 90]); axis tight; shading interp; colorbar; drawnow;
figure(3);
PlotField2D(N,x,y,Div2D(uproj,vproj)); view([0 90]); axis tight; shading interp; colorbar; drawnow;

figure(4);
PlotField2D(N,x,y,urot); view([0 90]); axis tight; shading interp; colorbar; drawnow;
figure(5);
PlotField2D(N,x,y,vrot); view([0 90]); axis tight; shading interp; colorbar; drawnow;

figure(6);
PlotField2D(N,x,y,Div2D(urot,vrot)); view([0 90]); axis tight; shading interp; colorbar; drawnow;

figure(7);
PlotField2D(N,x,y,Div2D(u,v)); view([0 90]); axis tight; shading interp; colorbar; drawnow;