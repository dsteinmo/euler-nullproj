%
% Drives the SMPM incompressible Euler code.  Sets up and solves a falling
% density blob problem.
%
% 7 July 2015.
% Sumedh Joshi
% Cornell University.

   % Set some constants / inputs.
   n    = 6;
   mx   = 5;
   mz   = 5;
   Lx   = 1.0;
   Lz   = 1.0;
   rho0 = 1000.0;
   g    = 9.8;

   % Set the type of pressure projection.
%   ptype = 'nullspace';
%   ptype = 'poisson';
   ptype = 'postproject';
%  ptype = 'none';

   % Build the operator matrices.
   r = n * n * mx * mz;
   [Dx Dz E0 E1 B0 B1] = smpm_assemble_2D_cartesian( n, mx, mz, Lx, Lz );

   % Get a penalty coefficient.

   % Build the Poisson matrix and its null space.
   EL = 1000 * (E0 + E1) + 60 * B1; % Need to set a penalty parameter here.
   L = Dx * Dx + Dz * Dz + EL;
   [u0, junk, junk] = svds( L, 1, 0 );

   % Build a cartesian grid.
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
   fprintf( 'Computing null space of the divergence operator.\n' );
   N = compute_divergence_nullspace( D, n, mx, mz );

   % Build the vector C0 continuity operator.
   E_C0 = [ (E0 + 0*E1 + B0) zeros(r,r); zeros(r,r) (E0 + 0*E1 + B0) ];

   % Set the regularization coefficient.
   tau = 1.0e4;

   % Solve for a truncated SVD.
   fprintf( 'Computing the SVD of normal equations.\n' );
   [U S V] = setup_nullspace_projection( tau, N, E_C0 );

   % Build an initial velocity.
   x0    = Lx/2;
   z0    = Lz/2;
   a     = Lx / 8;
   Umax  = 0.1;
   %ux0 = U * exp( -( ( x - x0 ).^2 + ( z - z0 ).^2 ) / a.^2 );
   ux0 = 0 * x;
   uz0 = 0 * ux0;

   % Build an initial density.
   rhoi = rho0 * exp( -( ( x - x0).^2 + (z - z0).^2 ) / (a/2).^2 );

   % Set some time-stepping parameters.
   t_final = 0.2;
   min_dt  = 2e-4;
   min_dx  = z(2) - z(1);  % XXX: Only works for cartesian grids.
   c       = sqrt(max( ux0.^2 + uz0.^2 ));
   dt      = min( min_dx / c / 10, min_dt );
   t       = linspace( 0, t_final, ceil( t_final / dt ) );

   % Allocate arrays to store the field variables.
   ux   = zeros( r, length(t) );
   uz   = zeros( r, length(t) );
   rho  = zeros( r, length(t) );
   Drho = zeros( r, length(t) ); % Estimate of material derivative of density.
   mass = zeros( length(t), 1 );    % Total mass.
   mass(1) = integrate_grid_function( rhoi, x, z, n, mx, mz );

   %
   % Begin time-stepping.
   ux(  :, 1 ) = ux0;
   uz(  :, 1 ) = uz0;
   rho( :, 1 ) = rhoi;
   for ii = 2:length(t)

      % Apply the non-linear advection operator.
      Aux = apply_smpm_advection( ux(:,ii-1), Dx, Dz, ux(:,ii-1), uz(:,ii-1), Lx, Lz, n, mx, mz );
      Auz = apply_smpm_advection( uz(:,ii-1), Dx, Dz, ux(:,ii-1), uz(:,ii-1), Lx, Lz, n, mx, mz );

      % Apply the advective operator to the density.
      Arho = apply_smpm_advection( rho(:,ii-1), Dx, Dz, ux(:,ii-1), uz(:,ii-1), Lx, Lz, n, mx, mz );

      % Update the current velocity.
      ux(:,ii) = ux(:,ii-1) + dt * Aux;
      uz(:,ii) = uz(:,ii-1) + dt * Auz - dt * g * rho(:,ii-1) / rho0 ;

      % Update the current density.
      rho(:,ii) = rho(:,ii-1) + dt * Arho;

      % Get an estimate of the material derivative of the density (for estimating satisfaction of the divergence-free-ness).
      Drho(:,ii) = ( rho(:,ii) - rho(:,ii-1) ) / dt + ( Dx * rho(:,ii) ) .* ux(:,ii) + ( Dz * rho(:,ii) ) .* uz(:,ii);

      % Estimate total mass.
      mass(ii) = integrate_grid_function( rho(:,ii), x, z, n, mx, mz );


      % Project onto the divergence-free basis.
      switch ptype
         case 'nullspace'
            [ ux(:,ii), uz(:,ii) ] = apply_nullspace_projection( ux(:,ii), uz(:,ii), N, U, S, V );
         case 'poisson'

            % Set up a right-hand-side.
            b = -D * [ ux(:,ii); uz(:,ii) ];
            b = b - u0 * u0' * b;

            % Solve the Poisson equation.
            p = L \ b;

            % Update the current velocities.
            Gp = G*p;
            ux(:,ii) = ux(:,ii) - dt * Gp(1:r);
            uz(:,ii) = uz(:,ii) - dt * Gp(r+1:end);

         case 'postproject'

            % Set up a right-hand-side.
            b = -D * [ ux(:,ii); uz(:,ii) ];
            b = b - u0 * u0' * b;

            % Solve the Poisson equation.
            p = L \ b;

            % Update the current velocities.
            Gp = G*p;
            ux(:,ii) = ux(:,ii) - dt * Gp(1:r);
            uz(:,ii) = uz(:,ii) - dt * Gp(r+1:end);

            % Project onto the divergence-free basis.
            iiu = [ux(:,ii); uz(:,ii)];
            iiu = N * ( N' * iiu );
            ux(:,ii) = iiu(1:r);
            uz(:,ii) = iiu(r+1:end);
      end

      % Print out a CFL number.
      iiCFL = sqrt( max( ux(:,ii).^2 + uz(:,ii).^2 ) ) * dt / min_dx;
      fprintf( [ 'Time step ' num2str( ii ) ', CFL number:', num2str( iiCFL ) '\n'] );

   end

   % Reshape some variables for plotting.
   X  = reshape( x, n * mx, n * mz );
   Z  = reshape( z, n * mx, n * mz );
   Ux = reshape( ux(:,end), n * mx, n * mz );
   Uz = reshape( uz(:,end), n * mx, n * mz );
   R  = reshape( rho(:,end), n * mx, n * mz );
