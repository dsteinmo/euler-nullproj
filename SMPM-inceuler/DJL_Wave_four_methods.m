%
% Drives the SMPM incompressible Euler code.
% Solves a sample DJL wave propagating inviscidly.
%
% 7 July 2015.
% Sumedh Joshi
% Cornell University.

   % Set some constants / inputs.
   n    = 5;
   mx   = 20;
   mz   = 10;
   Lx   = 10.0;
   Lz   = 0.15;
   rho0 = 1000.0;
   g    = 9.8;

   % Build the operator matrices.
   r = n * n * mx * mz;
   fprintf(['Assembling operator matrices.\n']);
   [Dx Dz E0 E1x E1z B0x B0z B1] = smpm_assemble_2D_cartesian( n, mx, mz, Lx, Lz );

   % Load the DJL wave data.
   load samplewave;

   % Set the bounds to be [0,Lz].
   zz = zz + Lz;

   % Get a penalty coefficient.

   % Build a cartesian grid.
   fprintf(['Building the grid.\n']);
   [x, z] = smpm_build_cartesian_mesh( n, mx, mz, [0 Lx], [0 Lz] );

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

   % Interpolate the DJL data onto the element grid we just built.
   ux0  = interp2( xx, zz, u, x, z,   'spline', 0.0 );
   uz0  = interp2( xx, zz, w, x, z,   'spline', 0.0 );
   rhoi = interp2( xx, zz, eta, x, z, 'spline', 0.0 );
   rhob = interp1( zz(:,1), rho(:,1), z, 'spline' );

      % Set the top and bottom of the domain.
      ux0 = reshape( ux0, n * mz, n * mx );
      uz0 = reshape( uz0, n * mz, n * mx );
      ux0( 1, : )   = ux0( 2, : );
      ux0( end, : ) = ux0( end-1, : );
      %uz0( 1, : )   = uz0( 2, : );
      %uz0( end, : ) = uz0( end-1, : );
      ux0 = ux0(:);
      uz0 = uz0(:);

   % Fix the units on the DJL data.
   rhob = rhob * 1000;
   rhoi = rhoi * 1000;

   % Build the vector C0 continuity operator.
   E_C0 = [ (E0 + E1x + E1z + B0x) zeros(r,r); zeros(r,r) (E0 + E1x + E1z + B0z) ];

   % Set the regularization/penalty coefficient.
   tau(1) = 1.0e4;
   tau(2) = 1.0e4;
   tau(3) = 1.0e4;
   tau(4) = 1.0e4;

   % Set some time-stepping parameters.
   t_final = 1.00;
   min_dt  = 0.05;
   min_dx  = z(2) - z(1);  % XXX: Only works for cartesian grids.
   c       = sqrt(max( ux0.^2 + uz0.^2 ));
   dt      = min( min_dx / c / 10, min_dt );

   % Solve the incompressible Euler equations three times, each with a different projection method.
   %ptypes = { 'postnull', 'poisson', 'postproject', 'nullspace-direct' };
   %fname  = { 'raytay_postnull', 'raytay_poisson', 'raytay_postproject', 'raytay_nullspace' };
   ptypes = { 'poisson' };
   fname  = { 'djl_poisson' };
   for ii = 1:length(ptypes)
      [ux uz rho t] = solve_incompressible_euler( n, mx, mz, x, z, ux0, uz0, rhoi, dt, t_final, ptypes{ii}, tau(ii) );
      save( fname{ii} );
   end
