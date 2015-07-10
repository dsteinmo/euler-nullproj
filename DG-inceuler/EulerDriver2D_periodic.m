%now trying hard cut-off filter, cutoff=.45
%and trying to just put in w=0 at boundaries initially.
%none of that really did anything.

%if this doesn't work, try no filter on rho next.
%=>no filter on rho helps a lot.
%removing initial field filters, too.

clear;
%close all;
TestAll_belize;
addpath /RAID/dsteinmo/localdivproj/
addpath /RAID/dsteinmo/bandLimFourierInterp/
addpath /RAID/dsteinmo/NUDG/Codes2Dquad
% Driver script for solving the 2D Euler equations
Globals2D;
vmapW = 0*vmapW;
vmapB = 0*vmapB;
vmapM = 0*vmapM;
vmapP = 0*vmapP;

%keyboard;
% Order of polynomials used for approximation 
N = 8; 


BCSolution = @NHSWMBC2D;


% read mesh from file
%[Nv, VX, VY, K, EToV, BCType,node,edge] = readmsh('../gups_xzplane2.msh');
[Nv, VX, VY, K, EToV, BCType,node,edge] = readmsh('../gups_xzplane3.msh');
%[Nv, VX, VY, K, EToV, BCType,node,edge] = readmsh('../../../gups_xzplane_coarse.msh');
%re-scale to unit square
VY=-VY;

g=9.81;

% set up nodes and basic operations
StartUp2Dquad;

% turn Neuman into walls
ids = find(BCType==7); 
BCType(ids) = Wall;

%need to strip off BC tags that are bad., otherwise poisson solver
%won't work:

BuildBCMaps2D;
vmapWold = vmapW;
mapWold = mapW;
%make x-direction periodic
%                      xmax-xmin, dummy value for y.
NODETOL=1e-12;
BuildPeriodicMaps2D(max(x(:))-min(x(:)),100);
vmapW = vmapB;
mapW = mapB;

vmapPer = setdiff(vmapWold,vmapW);
mapPer = setdiff(mapWold,mapW);

BCType = CorrectBCTablequad(EToV,find(abs(VY-min(y(:)))<NODETOL | abs(VY-max(y(:)))<NODETOL),Wall,K);
%figure(27);
%PlotDomain2D;
%return;

%Now just need to find the 4 corner points and set
%vmapM = vmapP in those. Do we need this?
if 0 ==1
corvmapM = find(x == min(x(:)) & y == min(y(:)));
corind = find(vmapM==corvmapM);
vmapP(corind) = vmapM(corind);

corvmapM = find(x == max(x(:)) & y == max(y(:)));
corind = find(vmapM == corvmapM);
vmapP(corind) = vmapM(corind);

corvmapM = find(x == max(x(:)) & y == min(y(:)));
corind = find(vmapM == corvmapM);
vmapP(corind) = vmapM(corind);

corvmapM = find(x == min(x(:)) & y == max(y(:)));
corind = find(vmapM == corvmapM);
vmapP(corind) = vmapM(corind);

%figure(16);
%subplot(2,1,1);
%plot(x(corvmapM),y(corvmapM),'.');
%subplot(2,1,2);
%plot(x(vmapP(corind)),y(vmapP(corind)),'.');
%return;

corind= 0*corind;
%faces are a bit weirder, more than 1 face can have the corner point.

cormapM = find(Fx == min(Fx(:)) & Fy == min(Fy(:)));
for j=1:length(cormapM);
    corind = find(mapM == cormapM(j));
    mapP(corind) = mapM(corind);
end

cormapM = find(Fx == max(Fx(:)) & Fy == max(Fy(:)));
for j=1:length(cormapM);
    corind = find(mapM == cormapM(j));
    mapP(corind) = mapM(corind);
end

cormapM = find(Fx == max(Fx(:)) & Fy == min(Fy(:)));
for j=1:length(cormapM);
    corind = find(mapM == cormapM(j));
    mapP(corind) = mapM(corind);
end

cormapM = find(Fx == min(Fx(:)) & Fy == max(Fy(:)));
for j=1:length(cormapM)
    corind = find(mapM == cormapM(j));
    mapP(corind) = mapM(corind);
end

end %endif 0==1

%figure(19); clf;
%PlotDomain2D;
%figure(20); clf;
%subplot(4,1,1);
%plot(x(vmapW),y(vmapW),'.');
%subplot(4,1,2);
%plot(Fx(mapW),Fy(mapW),'.');
%subplot(4,1,3);
%plot(x(vmapPer),y(vmapPer),'.');
%subplot(4,1,4);
%plot(Fx(mapPer),Fy(mapPer),'.');
%return;

%figure(21); clf;
%plot(x(corind),y(corind),'.');
%return;
close all;

% compute initial condition
Q = zeros(Np,K,4);

%load GUPS data
%load('../../../../../../samplewave.mat')
%clear eta %c

%define background density
%zpyc = 0.10;
%dpyc = 0.02;
%drho_halved = 0.02;
%zza=1/dpyc; zzb=zpyc/dpyc;
%rhobarDG = -drho_halved*tanh(zza*y-zzb);
%rhobar = -drho_halved*tanh(zza*zz+.05/dpyc);

%figure(1);
%subplot(3,1,1);
%pcolor(xx,zz,rho-1); shading flat; colorbar;
%subplot(3,1,2);
%pcolor(xx,zz,rhobar); shading flat; colorbar;
%subplot(3,1,3);
%pcolor(xx,zz,rho-1-rhobar); shading flat; colorbar;
%return;

%figure(1);
%pcolor(xx,zz,rho); shading flat; colorbar;  drawnow;
%sample-wave is hideously-well resolved, so downsampling
%for our grids with < 100,000 points is fine.

if 0==1

skip=4;

zz=-zz;
zz = flipud(zz);
%keyboard;
rhoprime = rho-1-rhobar;
%tmp = [-flipud(rhoprime) -flipud(fliplr(rhoprime)); rhoprime fliplr(rhoprime)];
tmp = [-flipud(rhoprime); rhoprime];
%return to normal by removing middle strip
xxe = [xx;
       xx];
%xxe = [xx  xx+max(xx(:));
%       xx  xx+max(xx(:))];
%zze = [zz zz; zz+min(zz(:)) zz+min(zz(:))]; %for negative z-coord
%zze = [zz+max(zz(:)) zz+max(zz(:)); zz zz];  %positive z-coord
zzup = zz+max(zz(:));
zze = [zzup;
       zz];
%figure(18);
%plot(xxe,zze,'.'); drawnow; return;

%could do it with fft's?
%tmphat = tmphat(end/4+1:3*end/4,end/4+1:3*end/4);
%tmp = real(ifft(tmphat));
%sz = size(tmp);
%tmp = (1024*512/sz(1)/sz(2));

%tmp = tmp(1:skip:end,1:skip:end);
%xxe = xxe(1:skip:end,1:skip:end);
%zze = zze(1:skip:end,1:skip:end);

%set-up low-pass filter
Nx = length(xxe(1,:));
Nz = length(zze(:,1));
Lx = max(xxe(:));
Lz = max(zze(:));
dk = 2*pi/Lx;  %Nyquist wavenumbers
dl = 2*pi/Lz;

k=[0:Nx/2-1 Nx/2 -Nx/2+1:-1]'*dk;
l=[0:Nz/2-1 Nz/2 -Nz/2+1:-1]'*dl;
[kk,ll]=meshgrid(k,l);

%Build filter
kmax=max(kk(:));
lmax=max(ll(:));
%% SUBICH
cutoff=0.45; %% 0.65 typically,  CUTOFF OF 0 CORRESPONDS TO HYPERVISCOSITY
kcrit=kmax*cutoff;
lcrit=lmax*cutoff;
f_order=2; % filter order
epsf=1e-15; %% FILTER STRENGTH AT HIGH WAVENUMBERS

myfilt=ones(size(kk));

%filter in x
mymask = (abs(kk)<kcrit);
%myfilt = myfilt.*(mymask + ... % Less than cutoff, no filtering
%           (1-mymask).* ... % exponential filter
%           exp(log(epsf)*((abs(kk)-kcrit)./(max(kk(:))-kcrit)).^f_order));
myfilt = myfilt.*(mymask + (1-mymask).*(kk*0));

%filter in y
mymask = (abs(ll)<lcrit);
%myfilt = myfilt.*(mymask + (1-mymask).* ...
%           exp(log(epsf).*((abs(ll)-lcrit)./(max(ll(:))-lcrit)).^f_order));
myfilt = myfilt.*(mymask + (1-mymask).*(ll*0));

% remove the Nyquist frequency entirely
myfilt(:,floor(Nx/2+1)) = 0;
myfilt(floor(Nz/2+1),:) = 0;

%figure(1);
%pcolor(kk,ll,myfilt); shading flat; drawnow; colorbar;
%return;

clear kk ll

%apply low-pass filter
figure(11);
%subplot(2,1,1);
pcolor(xxe,zze,tmp); shading flat; colorbar;
%tmp = real(ifft2(myfilt.*fft2(tmp)));
%subplot(2,1,2);
%pcolor(xxe,zze,tmp); shading flat; colorbar; drawnow;


%try band-limited interpolation
disp('interpolating rho to DG grid...');
rhoprimeDG = bandLimFourierInterp2D(xxe,zze,tmp,x(:),y(:),1e10);
rhoprimeDG = reshape(rhoprimeDG,Np,K);

figure(12);clf;
subplot(2,2,1);
pf2dquad(N,x,y,rhoprimeDG); shading interp; colorbar; drawnow;
title('\rho prime DG');
rhoDG = rhoprimeDG+rhobarDG;
subplot(2,2,2);
pf2dquad(N,x,y,rhoDG); shading interp; colorbar; drawnow;
title('\rho DG');
subplot(2,2,3);
pcolor(xx,zz,rhoprime); shading interp; colorbar; drawnow;
title('\rho prime spectral');
subplot(2,2,4);
pcolor(xx,zz,rho-1); shading interp; colorbar; drawnow;
title('\rho spectral');
%return;

tmp = [-flipud(eta); eta];

%tmp = [flipud(u); u];
%tmp = [flipud(u) fliplr(flipud(u)); u fliplr(u)]; %haven't tried this (even u)

%tmp = tmp(1:skip:end,1:skip:end);
figure(2); clf;
subplot(2,1,1);
pcolor(xxe,zze,tmp); shading flat; colorbar;
%tmp = real(ifft2(myfilt.*fft2(tmp)));

disp('interpolating eta to DG grid...');
etaDG = bandLimFourierInterp2D(xxe,zze,tmp,x(:),y(:),1e10);
etaDG = reshape(etaDG,Np,K);
%clear u;
subplot(2,1,2);
pf2dquad(N,x,y,etaDG); shading flat; colorbar; drawnow;
[etaxDG,etayDG] = Grad2D(etaDG);
uDG = c*etayDG;
wDG =-c*etaxDG; 

end %end if.

%close all;
%figure(1);
%pf2d(N,x,y,uDG); colorbar;
%drawnow;
%figure(2);
%pf2d(N,x,y,wDG); colorbar;
%drawnow;

%pause;

% figure(2);
% PlotField2D(N,x,y,rhoDG); view([0 90]); axis tight; shading flat; colorbar; drawnow; return;
% return;


%save('../gups_DG_N4_modinterp.mat','rhoDG','uDG','wDG');
load('../gups_DG_N8_modinterp.mat');

%uDG(uDG<1e-6) =0;
%wDG(wDG<1e-7) =0;

rho = rhoDG;
u = uDG;
v = wDG;

Lx =max(x(:))-min(x(:));
Ly = max(y(:))-min(y(:));

%u = u.*sin(pi*x/Lx);
%v = v.*sin(pi*y/Ly);

%figure(2);
%pf2d(N,x,y,rho); shading interp; colorbar;

clear rhoDG uDG wDG


%rho = 0.003*(-tanh(50*(y+.005*x-.1))+1)/2;
%figure(1); clf;
%pf2dquad(N,x,y,rho); colorbar; drawnow; return;
%Q(:,:,1) = rho; %density perturbation, not full density.
%Q(:,:,2) = 0*rho;
%Q(:,:,3) = 0*rho;
%Q(:,:,4) = 0*rho;

%a0=0.02;
%lambda=.35;
%eta = a0*sech((x-0.5*Lx)/lambda).^2.*sin(pi*y/Ly);
%rhobar = 0.003*(-tanh(50*(y-.1))+1)/2;
%rho    = 0.003*(-tanh(50*(y-.1+eta))+1)/2;

%figure(1); clf;
%subplot(3,1,1);
%pf2dquad(N,x,y,rhobar); colorbar;
%subplot(3,1,2);
%pf2dquad(N,x,y,eta); colorbar;
%subplot(3,1,3);
%pf2dquad(N,x,y,rho); colorbar;
%return;
%c0=0.05; %wave speed
%u =-(pi*a0*c0/Ly)*sech((x-0.5*Lx)/lambda).^2.*cos(pi*y/Ly);
%w = (2*a0*c/lambda)*sech((x-0.5*Lx)/lambda).^2.*tanh((x-.5*Lx)/lambda).*sin(pi*y/Ly);;

%figure(1); clf;
%subplot(3,1,1);
%pf2dquad(N,x,y,rho); colorbar;
%subplot(3,1,2);
%pf2dquad(N,x,y,u); colorbar;
%subplot(3,1,3);
%pf2dquad(N,x,y,w); colorbar;
%return;

%try this to force the DG grid's bc's:
%v(vmapW)=0;

%spacefilt=.5*(tanh(6*(x-0.25*Lx))+1) - .5*(tanh(6*(x-0.75*Lx))+1);
%u = u.*spacefilt;
%v = v.*spacefilt;
%v = v.*(tanh(x-0.5*Lx)+1)/2;

%rho = rho-1;

Q(:,:,1) = rho;
Q(:,:,2) = u;
Q(:,:,3) = v;
Q(:,:,4) = 0*rho;



%figure(1); clf;
%pf2dquad(N,x,y,u); shading interp; colorbar; axis tight;
%drawnow;
%keyboard
close all;
% Solve Problem
FinalTime = 1000; %48 = roughly tank traversal time
%[Q] = Euler2Dpressure(Q, FinalTime, BCSolution,g); 
[Q] = Euler2Dpresstiming(Q,2,BCSolution,g); %for timing test.
