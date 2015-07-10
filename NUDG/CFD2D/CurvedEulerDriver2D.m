% Driver script for solving the 2D vacuum Euler's equations 
Globals2D;

% Order of polynomial approximation (N) 
N = 8;

% set up simulation type
sim = 'IsentropicVortex'; 
%sim = 'ChannelFlow';

%fluxtype = 'Roe';
%fluxtype = 'LF';
fluxtype = 'HLL';

switch sim
case {'IsentropicVortex'}
  filename = 'vortexA04.neu';
  InitialSolution = @IsentropicVortexIC2D;
  ExactSolution   = @IsentropicVortexIC2D;
  BCSolution      = @IsentropicVortexBC2D;
case {'ChannelFlow'}
  filename = 'Euler01.neu';
  InitialSolution = @ChannelIC2D;
  ExactSolution   = [];
  BCSolution      = @ChannelBC2D;
otherwise 
  disp('Simulation case unknown');  stop;
end

filename = 'CNScylK930.neu';
% Read in Mesh
[Nv, VX, VY, K, EToV, BCType] = MeshReaderGambitBC2D(filename);

% Initialize solver and construct grid and metric
StartUp2D;

% adjust curved elements according to simulation type
% hard wired, cylinder (centered at (0,0) with radius .5)
[k,f] = find(BCType==Cyl);
curved = [];
if(~isempty(k))
 cylfaces = [k,f];
 curved = sort(unique(k)); 
 MakeCylinder2D_derek(cylfaces, 1, 0, 0);
 % turn cylinders into walls
 ids = find(BCType==Cyl);  BCType(ids) = Wall;
end
straight = setdiff(1:K, curved);
BuildBCMaps2D

% compute initial condition (time=0)
Q = feval(InitialSolution, x, y, 0);

% Solve Problem
FinalTime = 10;
[Q] = CurvedEuler2D(Q, FinalTime, ExactSolution, BCSolution, fluxtype); 
