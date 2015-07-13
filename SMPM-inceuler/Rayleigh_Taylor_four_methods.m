%
% Drives the SMPM incompressible Euler code.  Sets up and solves a falling
% density blob problem that is something like a Rayleigh-Taylor instability.
%
% 7 July 2015.
% Sumedh Joshi
% Cornell University.

   % Set some constants / inputs.
   n    = 8;
   mx   = 10;
   mz   = 10;
   Lx   = 1.0;
   Lz   = 1.0;
   rho0 = 1000.0;
   g    = 9.8;

   % Build the operator matrices.
   r = n * n * mx * mz;
   fprintf(['Assembling operator matrices.\n']);
   [Dx Dz E0 E1x E1z B0x B0z B1] = smpm_assemble_2D_cartesian( n, mx, mz, Lx, Lz );

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
   %N = compute_divergence_nullspace( D, n, mx, mz );

   % Build the vector C0 continuity operator.
   E_C0 = [ (E0 + E1x + E1z + B0x) zeros(r,r); zeros(r,r) (E0 + E1x + E1z + B0z) ];

   % Set the regularization/penalty coefficient.
   tau(1) = 1.0e4;
   tau(2) = 1.0e4;
   tau(3) = 1.0e4;
   tau(4) = 1.0e4;

   % Build an initial velocity.
   x0    = Lx/2;
   z0    = Lz/2;
   a     = Lx / 8;
   Umax  = 0.1;
   ux0 = 0 * x;
   uz0 = 0 * ux0;

   % Build an initial density.
   rhoi = ( 0.1 * rho0 ) * exp( -( ( x - x0).^2 + (z - z0).^2 ) / (a/2).^2 );

   % Set some time-stepping parameters.
   t_final = 1.00;
   min_dt  = 0.05;
   min_dx  = z(2) - z(1);  % XXX: Only works for cartesian grids.
   c       = sqrt(max( ux0.^2 + uz0.^2 ));
   dt      = min( min_dx / c / 10, min_dt );

   % Solve the incompressible Euler equations three times, each with a different projection method.
   ptypes = { 'postnull', 'poisson', 'postproject', 'nullspace-direct' };
   fname  = { 'raytay_postnull', 'raytay_poisson', 'raytay_postproject', 'raytay_nullspace' };
   for ii = 1:length(ptypes)
      [ux uz rho t] = solve_incompressible_euler( n, mx, mz, x, z, ux0, uz0, rhoi, dt, t_final, ptypes{ii}, tau(ii) );
      save( fname{ii} );
   end
