%
% Drives the SMPM incompressible Euler code.  Sets up and solves a falling
% density blob problem.
%
% 7 July 2015.
% Sumedh Joshi
% Cornell University.

   % Set some constants / inputs.
   n    = 10;
   mx   = 10;
   mz   = 10;
   Lx   = 1.0;
   Lz   = 1.0;
   rho0 = 1000.0;
   g    = 9.8;

   % Set the type of pressure projection.
%   ptype = 'nullspace-direct';
   ptype = 'nullspace-iterative';
%  ptype = 'poisson';
%   ptype = 'postproject';
%   ptype = 'none';

   % Build the operator matrices.
   r = n * n * mx * mz;
   fprintf(['Assembling operator matrices.\n']);
   [Dx Dz E0 E1 B0 B1] = smpm_assemble_2D_cartesian( n, mx, mz, Lx, Lz );

   % Get a penalty coefficient.

   % Build a cartesian grid.
   fprintf(['Building the grid.\n']);
   [x, z] = smpm_build_cartesian_mesh( n, mx, mz, [0 1], [0 1] );

   % Build the divergence and gradient operators.

      % Divergence.
      r  = n * n * mx * mz;
      D = sparse( r, 2 * r );
      D(:,1:r)     = Dx;
      D(:,r+1:end) = Dz;

      % Gradient.
      G = sparse( 2 * r, r );
      G( 1:r, : )     = Dx;
      G( r+1:end, : ) = Dz;

   % Assemble the null space of the divergence operator.
   %N = null( full( D ) ); % This is slow and dumb.
   N = compute_divergence_nullspace( D, n, mx, mz );

   % Build the vector C0 continuity operator.
   E_C0 = [ (E0 + 0*E1 + B0) zeros(r,r); zeros(r,r) (E0 + 0*E1 + B0) ];

   % Set the regularization coefficient.
   tau = 1.0e4;

   % Build an initial velocity.
   x0    = Lx/2;
   z0    = Lz/2;
   a     = Lx / 8;
   Umax  = 0.1;
   ux0 = 0 * x;
   uz0 = 0 * ux0;

   % Build an initial density.
   rhoi = rho0 * exp( -( ( x - x0).^2 + (z - z0).^2 ) / (a/2).^2 );

   % Set some time-stepping parameters.
   t_final = 0.25;
   min_dt  = 1e-3;
   min_dx  = z(2) - z(1);  % XXX: Only works for cartesian grids.
   c       = sqrt(max( ux0.^2 + uz0.^2 ));
   dt      = min( min_dx / c / 10, min_dt );
   t       = linspace( 0, t_final, ceil( t_final / dt ) );

   % Solve the incompressible Euler equations.
   [ux uz rho] = solve_incompressible_euler( n, mx, mz, x, z, ux0, uz0, rhoi, dt, t_final, ptype, tau );
