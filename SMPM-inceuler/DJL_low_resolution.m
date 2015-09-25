%
% Drives the SMPM incompressible Euler code.
% Solves a sample DJL wave propagating inviscidly.
%
% 7 July 2015.
% Sumedh Joshi
% Cornell University.

   % Set some constants / inputs.
   n    = 10;
   mx   = 8;
   mz   = 3;
   Lx   = 12.0;
   Lz   = 0.15;
   rho0 = 1000.0;
   g    = 9.8;

   % Set some time-stepping parameters.
   t_final = 60.0;
   min_dt  = 0.05;

   % Set the types of simulations to run.
   ptypes = { 'poisson', 'postproject', 'postnormal', 'postnull' };
   fname  = { 'djl_poisson', 'djl_postproject', 'djl_postnormal', 'djl_postnull' };

   % Set the regularization/penalty coefficient.
   tau(1) = 1.0e4;
   tau(2) = 1.0e4;
   tau(3) = 1.0e4;
   tau(4) = 1.0e4;

   % Data reduction factor.
   reduce = 100;

   %% Below here there are no inputs needing to be changed.

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
   rhoi = interp2( xx, zz, rho - repmat(rho(:,1), 1, size(rho,2)), x, z, 'spline', 0.0 );
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
   rhob = rhob * 1000 - rho0;
   rhoi = rhoi * 1000;

   % Build the vector C0 continuity operator.
   E_C0 = [ (E0 + 0*E1x + 0*E1z + B0x) zeros(r,r); zeros(r,r) (E0 + 0*E1x + 0*E1z + B0z) ];

   % Set some time-stepping parameters.
   min_dx  = min( z(2) - z(1), x(n*mz + 1) - x(1) );  % XXX: Only works for cartesian grids.
   c       = sqrt(max( ux0.^2 + uz0.^2 )) + 0 * c;
   dt      = min_dt;

   % Solve the incompressible Euler equations three times, each with a different projection method.
   for ii = 1:length(ptypes)
         [ux uz rho t] = solve_incompressible_euler( n, mx, mz, x, z, ux0, uz0, rhoi, rhob, dt, t_final, ptypes{ii}, tau(ii) );

         % Reshape the outputs.
         rho = reshape( rho, [ n * mz, n * mx, length(t) ] );
         ux  = reshape( ux, [ n * mz, n * mx, length(t) ] );
         uz  = reshape( uz, [ n * mz, n * mx, length(t) ] );

         % Decimate the data prior to saving it.
         if length(t) > reduce
            rho = rho( :, :, 1:reduce:end );
            ux  =  ux( :, :, 1:reduce:end );
            uz  =  uz( :, :, 1:reduce:end );
            try
               uxu = uxu( :, :, 1:reduce:end );
               uzu = uzu( :, :, 1:reduce:end );
            end
            t   =   t( 1:reduce:end );
         end

         % Save data.
         save( fname{ii} );
   end

   % Reshape ancillary arrays.
   rhob = reshape( rhob, n * mz, n * mx );
   x    = reshape( x, n * mz, n * mx );
   z    = reshape( z, n * mz, n * mx );

   % Check some energy stuff.
   E  = 0.5 * rho.*(ux.^2 + uz.^2);
   Eu = 0.5 * rho.*(uxu.^2 + uzu.^2);
   Etot  = integrate_grid_function( reshape( E,  n*n*mx*mz, length(t)), x(:), z(:), n, mx, mz );
   Eutot = integrate_grid_function( reshape( Eu, n*n*mx*mz, length(t)), x(:), z(:), n, mx, mz );

