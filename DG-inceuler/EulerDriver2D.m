%now try just averaging across the normal direction.
%==> averaging winning with t=28.5

%now trying to average both u and v
%==> worst one yet!

%now turning filtering back on, increasing dt by 10x. seeing what we get.

clear;
close all;

addpath "../"
setNUDGpaths;

% Driver script for solving the 2D Euler equations
Globals2D;
[Nv,VX, VY, K, EToV,BCType,node,edge] = readmsh('mysquare_refine.msh');

N=4;
% Read in Mesh


BCSolution = @NHSWMBC2D;


g=9.81;

% set up nodes and basic operations
StartUp2Dquad;

%h-refine
%done h-refine

% turn Neuman into walls
ids = find(BCType==7); 
BCType(ids) = Wall;

BuildBCMaps2D;

% compute initial condition
Q = zeros(Np,K,4);
rho = 1 + 0.003*(-tanh(5*(y+.1*x-.6))+1)/2 ; %this is good pycnocline for football
u = zeros(Np,K);
v = u;


%Old IC'S:
%rho = 1 + 0.03*(-tanh(20*(y-.5+.2*cos(2*pi*x)))+1)/2 + 1e-5*rand(Np,K);
%u = -(y-.5).*exp(-10*((x-.5).^2+(y-.5).^2));
%v = (x-.5).*exp(-10*((x-.5).^2+(y-.5).^2));

Q(:,:,1) = rho-1; %try just density perturbation
Q(:,:,2) = u;
Q(:,:,3) = v;
Q(:,:,4) = 0*rho;

FinalTime = 200;
[Q] = Euler2Dpressure(Q, FinalTime, BCSolution,g); %pressure projection implementation
