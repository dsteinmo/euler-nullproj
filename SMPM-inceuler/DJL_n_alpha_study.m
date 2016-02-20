%
% Drives the SMPM incompressible Euler code.
% Solves a sample DJL wave propagating inviscidly.
%
% 7 July 2015.
% Sumedh Joshi
% Cornell University.

   % Set some constants / inputs.
   n    = 12;
   mx   = 8;
   mz   = 3;
   Lx   = 12.0;
   Lz   = 0.15;
   rho0 = 1000.0;
   g    = 9.8;

   % Load the DJL wave data.
   load samplewave;
   c_wave = c;

   % Set some time-stepping parameters so that the wave propagates an integer number of element lengths.
   number_of_elements = 1;
   hx = Lx / mx;
   t_final = hx * number_of_elements / c_wave;
   min_dt  = 0.15;

   % Set the regularization/penalty coefficient array we want to study.
   tau    = 10^4;
   alpha1 = 10.^[1:5];
   alpha2 = 10.^[1:5];

   % Set the values of n we want to study.
   N = [4 6 8 10 12];

   % Data reduction factor.
   reduce = 100;

   %% Below here there are no inputs needing to be changed.

   % Build the operator matrices.
   r = n * n * mx * mz;
   fprintf(['Assembling operator matrices.\n']);
   [Dx Dz E0 E1x E1z B0x B0z B1] = smpm_assemble_2D_cartesian( n, mx, mz, Lx, Lz );

   % Set the bounds to be [0,Lz].
   zz = zz + Lz;

   % Get a penalty coefficient.


   % Loop over values of n.
   for kk = 1:length(N)

      % Set n.
      n = N(kk);

      % Build a cartesian grid.
      fprintf(['Building the grid.\n']);
      [x, z] = smpm_build_cartesian_mesh( n, mx, mz, [0 Lx], [0 Lz] );

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
      ux0 = ux0(:);
      uz0 = uz0(:);

      % Fix the units on the DJL data.
      rhob = rhob * 1000 - rho0;
      rhoi = rhoi * 1000;

      % Set some time-stepping parameters.
      min_dx  = min( z(2) - z(1), x(n*mz + 1) - x(1) );  % XXX: Only works for cartesian grids.
      c       = sqrt(max( ux0.^2 + uz0.^2 )) + 0 * c;
      dt      = min_dt;

         % Loop over alpha1's.
         for ii = 1:length(alpha1)

               % Loop over alpha2's.
               for jj = 1:length(alpha2)

                  % Solve the problem.
                  [ux uz rho t] = solve_incompressible_euler_wnp( n, mx, mz, x(:), z(:), ...
                                                                  ux0(:), uz0(:), rhoi(:), rhob(:), ...
                                                                  dt, t_final, ...
                                                                  tau, alpha1(ii), alpha2(jj) );

                  % Reshape the outputs.
                  rho = reshape( rho, [ n * mz, n * mx, length(t) ] );
                  ux  = reshape( ux, [ n * mz, n * mx, length(t) ] );
                  uz  = reshape( uz, [ n * mz, n * mx, length(t) ] );

                  % Get the final time.
                  rho = rho(:,:,end);
                  ux  = ux(:,:,end);
                  uz  = uz(:,:,end);

                  % Reshape ancillary arrays.
                  rhob = reshape( rhob, n * mz, n * mx );
                  x    = reshape( x,    n * mz, n * mx );
                  z    = reshape( z,    n * mz, n * mx );

                  % Reshape initial conditions.
                  ux0  = reshape( ux0, n * mz, n * mx );
                  uz0  = reshape( uz0, n * mz, n * mx );
                  rhoi = reshape( rhoi, n * mz, n * mx );

                  % Interpolate the initial conditions onto the last time.
                  ndx = number_of_elements * n;
                  diff_ux  = ux( :, ndx + 1:end ) - ux0( :, 1:end - ndx );
                  diff_uz  = uz( :, ndx + 1:end ) - uz0( :, 1:end - ndx );
                  diff_rho = rho( :, ndx + 1:end ) - rhoi( :, 1:end - ndx );

                  % Compare the solution at the final time to the initial time.
                  err_rel_u(ii,jj,kk)   = norm( diff_ux ) / norm( ux0(:, 1:end - ndx ) );
                  err_rel_w(ii,jj,kk)   = norm( diff_uz ) / norm( uz0(:, 1:end - ndx ) );
                  err_rel_rho(ii,jj,kk) = norm ( diff_rho ) / norm( rhoi(:, 1:end - ndx ) );;
                  err_u(ii,jj,kk)   = norm( diff_ux );
                  err_w(ii,jj,kk)   = norm( diff_uz );
                  err_rho(ii,jj,kk) = norm( diff_rho );

                  keyboard;

               end
         end
   end

