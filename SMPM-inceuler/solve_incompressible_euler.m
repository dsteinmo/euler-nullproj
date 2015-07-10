function [ux uz rho] = solve_incompressible_euler( n, mx, mz, x, z, ux0, uz0, rhoi, min_dt, t_final, ptype, tau )
% [ux uz rho] = solve_incompressible_euler( n, mx, mz, x, z, ux0, uz0, rho0, dt, t_final, ptype, tau );
%
%  Solves the incompressible Euler equations with a spectral multidomain
%  penalty method model on a cartesian (plaid) grid with an active density
%  and gravity.
%
%  Takes 12 arguments:
%
%     n, mx, mz      - grid discretization constants.
%     x, z           - grid generated by smpm_build_cartesian_grid.
%     ux0, uz0, rho0 - initial conditions, grid functions on x, z.
%     dt, t_final    - maximum timestep, final time.
%     ptype          - incompressibility enforcement method:
%                      {'poisson','nullspace-direct','nullspace-iterative','postproject','none'}.
%     tau            - either the regularization coefficient, or
%                      the factor multiplying the penalty coefficient.
%
%  Returns 3 arguments:
%
%     ux, uz, rho    - solution of the incompressible Euler equations.
%
%  7 July 2015
%  Sumedh Joshi
%  Cornell University

   % Set some constants / inputs.
   Lx   = max(x) - min(x);
   Lz   = max(z) - min(z);
   rho0 = 1000.0;
   g    = 9.8;

   % Build the operator matrices.
   r = n * n * mx * mz;
   [Dx Dz E0 E1 B0 B1] = smpm_assemble_2D_cartesian( n, mx, mz, Lx, Lz );

   % Get a penalty coefficient.

   % Build the Poisson matrix and its null space.
   if strcmp( ptype, 'poisson' ) || strcmp( ptype, 'postproject' )
      omega = 2.0 / ( n - 1 ) / n;
      kappa = omega;
      pen   = 1.0 / omega * ( 1.0 + ( 2 * kappa ) - ( 2 * sqrt( kappa^2 + kappa ) ) );
      EL = tau * pen * (E0 + E1) + pen^2 * B1;
      L = Dx * Dx + Dz * Dz + tau * EL;
      [u0, junk, junk] = svds( L, 1, 0 );
   end

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
   fprintf( 'Computing null space of the divergence operator.\n' );
   N = compute_divergence_nullspace( D, n, mx, mz );

   % Build the vector C0 continuity operator.
   E_C0 = [ (E0 + 0*E1 + B0) zeros(r,r); zeros(r,r) (E0 + 0*E1 + B0) ];

   % Solve for a truncated SVD if asked to do so.
   if strcmp( ptype, 'nullspace-direct' )
      fprintf( 'Computing the SVD of normal equations.\n' );
      [U S V] = setup_nullspace_projection( tau, N, E_C0 );
   end

   % Set up the normal equations for an interative method if asked to do so.
   if strcmp( ptype, 'nullspace-iterative' )
      fprintf( 'Assembling the normal equations.\n' );
      T = N'*N + ( E_C0 * N )' * ( E_C0 * N );
   end

   % Set some time-stepping parameters.
   min_dx  = z(2) - z(1);  % XXX: Only works for cartesian grids.
   c       = sqrt(max( ux0.^2 + uz0.^2 ));
   dt      = min( min_dx / c / 10, min_dt );
   t       = linspace( 0, t_final, ceil( t_final / dt ) );

   % Allocate arrays to store the field variables.
   ux   = zeros( r, length(t) );
   uz   = zeros( r, length(t) );
   rho  = zeros( r, length(t) );

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

      % Project onto the divergence-free basis.
      switch ptype

         case 'nullspace-direct'

            [ ux(:,ii), uz(:,ii) ] = apply_nullspace_projection( ux(:,ii), uz(:,ii), N, U, S, V );

         case 'nullspace-iterative'

            % Set up a right-hand-side.
            b = N' * [ ux(:,ii); uz(:,ii) ];

            % Solve with GMRES-Householder.
            [lambda err m] = compute_gmres_householder( T, b, b, 1e-9, r );

            % Construct the updated velocity.
            unew = N * lambda;
            ux(:,ii) = unew(1:r);
            uz(:,ii) = unew(r+1:end);

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
            Gp = G * p;
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

end
