
clear;
close all;
cd ../../
setNUDGpaths;
cd DG-inceuler/djl_wave
addpath('../../bandLimFourierInterp/')

% Driver script for solving the 2D Euler equations
Globals2D;
vmapW = 0*vmapW;
vmapB = 0*vmapB;
vmapM = 0*vmapM;
vmapP = 0*vmapP;

%keyboard;
% Order of polynomials used for approximation 
N = 4; 

%Do we want to interpolate from Fourier grid to DG grid, or use a cached file?
DO_INTERP = true;
BCSolution = @NHSWMBC2D;


% read mesh from file
[Nv, VX, VY, K, EToV, BCType,node,edge] = readmsh('gups_xzplane3.msh');
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
BuildPeriodicMaps2D(max(x(:))-min(x(:)),100);
vmapW = vmapB;
mapW = mapB;

vmapPer = setdiff(vmapWold,vmapW);
mapPer = setdiff(mapWold,mapW);

BCType = CorrectBCTablequad(EToV,find(abs(VY-min(y(:)))<NODETOL | abs(VY-max(y(:)))<NODETOL),Wall,K);
close all;

% compute initial condition
Q = zeros(Np,K,4);

%load GUPS data
load('samplewave.mat')
%clear eta %c

%define background density
zpyc = 0.10;
dpyc = 0.02;
drho_halved = 0.02;
zza=1/dpyc; zzb=zpyc/dpyc;
rhobarDG = -drho_halved*tanh(zza*y-zzb);
rhobar = -drho_halved*tanh(zza*zz+.05/dpyc);

if DO_INTERP==true

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

    save('gups_DG_N4_modinterp.mat','rhoDG','uDG','wDG');

end %end if.



load('gups_DG_N4_modinterp.mat');

rho = rhoDG;
u = uDG;
v = wDG;

Lx =max(x(:))-min(x(:));
Ly = max(y(:))-min(y(:));

clear rhoDG uDG wDG

Q(:,:,1) = rho;
Q(:,:,2) = u;
Q(:,:,3) = v;
Q(:,:,4) = 0*rho;



close all;
FinalTime = 1000; %48 = roughly tank traversal time
[Q] = Euler2Dpressure(Q, FinalTime, BCSolution,g); 

