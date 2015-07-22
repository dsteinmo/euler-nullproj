function [ux uz rho t] = solve_incompressible_euler( n, mx, mz, x, z, ux0, uz0, rhoi, rhob, max_dt, t_final, ptype, tau )
% [ux uz rho t] = solve_incompressible_euler( n, mx, mz, x, z, ux0, uz0, rho0, rhob, dt, t_final, ptype, tau );
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
%     rhob           - boussinesq density, or 0.0 if there is no background stratification.
%     dt, t_final    - maximum timestep, final time.
%     ptype          - incompressibility enforcement method:
%                      {'poisson','nullspace-direct','nullspace-iterative','postproject',postnull,'none'}.
%     tau            - either the regularization coefficient, or
%                      the factor multiplying the penalty coefficient.
%
%  Returns 3 arguments:
%
%     ux, uz, rho    - solution of the incompressible Euler equations.
%     t              - time.
%
%  7 July 2015
%  Sumedh Joshi
%  Cornell University

   % Set some constants / inputs.
   Lx      = max(x) - min(x);
   Lz      = max(z) - min(z);
   rho0    = 1000.0;
   g       = 9.8;
   CFL_MAX = 0.25;

   % Build the operator matrices.
   r = n * n * mx * mz;
   [Dx Dz E0 E1x E1z B0x B0z B1] = smpm_assemble_2D_cartesian( n, mx, mz, Lx, Lz );

   % Do some argument handling.
   if length( rhob ) == 1
      rhob = 0 * rhoi;
   end

   % Get a penalty coefficient.

   % Build the Poisson matrix and its null space.
   if strcmp( ptype, 'poisson' ) || strcmp( ptype, 'postproject' ) || strcmp( ptype, 'postnull' )
      omega = 2.0 / ( n - 1 ) / n;
      kappa = omega;
      pen   = 1.0 / omega * ( 1.0 + ( 2 * kappa ) - ( 2 * sqrt( kappa^2 + kappa ) ) );
      EL = tau * pen * (E0 + E1x + E1z) + pen^2 * B1;
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

   % Compute the derivative of the density stratification.
   rhob_z = Dz * rhob;

   % Assemble the null space of the divergence operator.
   fprintf( 'Computing null space of the divergence operator.\n' );
   N = compute_divergence_nullspace( D, n, mx, mz );

   % Build the vector C0 continuity operator.
   E_C0 = [ (E0 + E1x + E1z + B0x) zeros(r,r); zeros(r,r) (E0 + E1x + E1z + B0z) ];

   % Solve for a truncated SVD if asked to do so.
   if strcmp( ptype, 'nullspace-direct' ) || strcmp( ptype, 'postnull' )
      fprintf( 'Computing the SVD of normal equations.\n' );
      [U S V] = setup_nullspace_projection( tau, N, E_C0 );
   end

   % Set up the normal equations for an interative method if asked to do so.
   if strcmp( ptype, 'nullspace-iterative' )
      fprintf( 'Assembling the normal equations.\n' );
      T = N'*N + tau * ( E_C0 * N )' * ( E_C0 * N );
   end

   % Set some time-stepping parameters.
   min_dx  = z(2) - z(1);  % XXX: Only works for cartesian grids.
   c       = sqrt(max( ux0.^2 + uz0.^2 ));
   dt      = min( min_dx / c / 10, max_dt );
   t(1)    = 0.0;

   % Allocate arrays to store the field variables.
   ux   = zeros( r, 1 );
   uz   = zeros( r, 1 );
   rho  = zeros( r, 1 );

   %
   % Begin time-stepping.
   ux(  :, 1 ) = ux0;
   uz(  :, 1 ) = uz0;
   rho( :, 1 ) = rhoi;
   uzu = uz;
   uxu = ux;
   Gpx = 0 * ux(: ,1);
   Gpz = 0 * uz(:, 1);
   while t(end) < t_final

      % Apply the non-linear advection operator.
      Aux = apply_smpm_advection( ux(:,end), Dx, Dz, ux(:,end), uz(:,end), Lx, Lz, n, mx, mz );
      Auz = apply_smpm_advection( uz(:,end), Dx, Dz, ux(:,end), uz(:,end), Lx, Lz, n, mx, mz );

      % Apply the pressure gradient term.
      Aux = Aux - Gpx;
      Auz = Auz - Gpz;

      % Apply the advective operator to the density.
      Arho = apply_smpm_advection( rho(:,end), Dx, Dz, ux(:,end), uz(:,end), Lx, Lz, n, mx, mz );
      Arho = Arho - uz(:, end) .* rhob_z;

      % Update the current velocity.
      iiux = ux(:,end) + dt * Aux;
<<<<<<< HEAD
      iiuz = uz(:,end) + dt * Auz - dt * g * rho(:,end) ./ ( 0 * rhob + rho0 ) ;
=======
      iiuz = uz(:,end) + dt * Auz - dt * g * rho(:,end) ./ ( rho0 ) ;
>>>>>>> 0065346126cf12b7b6df5143f641e314392c6fd1

      % Update the current density.
      iirho = rho(:,end) + dt * Arho;

      uxu  = [uxu iiux];
      uzu  = [uzu iiuz];

      % Project onto the divergence-free basis.
      switch ptype

         case 'nullspace-direct'

            [ iiux, iiuz ] = apply_nullspace_projection( iiux, iiuz, N, U, S, V );

         case 'nullspace-iterative'

            % Set up a right-hand-side.
            b = N' * [ iiux; iiuz ];

            % Solve with GMRES-Householder.
            [lambda err m] = compute_gmres_householder( T, b, b, 1e-9, r );

            % Construct the updated velocity.
            unew = N * lambda;
            iiux = unew(1:r);
            iiuz = unew(r+1:end);

            % Display some statistics.
            fprintf([ '   GMRES converged in ', num2str( m ) ' iterations.\n'] );

         case 'poisson'

            % Set up a right-hand-side.
            b = -D * [ iiux; iiuz ] / dt;
            b = b - u0 * u0' * b;

            % Solve the Poisson equation.
            p = L \ b;

            % Update the current velocities.
            Gp = G*p;
            iiux = iiux + dt * Gp(1:r);
            iiuz = iiuz + dt * Gp(r+1:end);

         case 'postproject'

            % Set up a right-hand-side.
            b = -D * [ iiux; iiuz ] / dt;
            b = b - u0 * u0' * b;

            % Solve the Poisson equation.
            p = L \ b;

            % Update the current velocities.
            Gp = G * p;
            iiux = iiux + dt * Gp(1:r);
            iiuz = iiuz + dt * Gp(r+1:end);

            % Project onto the divergence-free basis.
            iiu = [iiux; iiuz];
            iiu = N * ( N' * iiu );
            iiux = iiu(1:r);
            iiuz = iiu(r+1:end);

         case 'postnull'

            % Set up a right-hand-side.
            b = -D * [ iiux; iiuz ] / dt;
            b = b - u0 * u0' * b;

            % Solve the Poisson equation.
            p = L \ b;

            % Update the current velocities.
            Gp = G * p;
            iiux = iiux + dt * Gp(1:r);
            iiuz = iiuz + dt * Gp(r+1:end);

            % Project onto the weakly continuous div-free basis.
            [ iiux, iiuz ] = apply_nullspace_projection( iiux, iiuz, N, U, S, V );

%            % Set up a right-hand-side.
%            b = N' * [ iiux; iiuz ];
%
%            % Solve with GMRES-Householder.
%            [lambda err m] = compute_gmres_householder( T, b, b, 1e-9, r );
%
%            % Construct the updated velocity.
%            unew = N * lambda;
%            iiux = unew(1:r);
%            iiuz = unew(r+1:end);
%
%            % Display some statistics.
%            fprintf([ '   GMRES converged in ', num2str( m ) ' iterations.\n'] );

      end

      % Print out a CFL number.
      umax  = sqrt( max( iiux.^2 + iiuz.^2 ) );
      iiCFL = sqrt( max( iiux.^2 + iiuz.^2 ) ) * dt / min_dx;
      fprintf( [ 'Time :' num2str(t(end)) ', CFL number:', num2str( iiCFL ) '\n'] );

      % Append the current time-step data.
      ux  = [ux iiux];
      uz  = [uz iiuz];
      rho = [rho iirho];

      % Set a new time-step.
      dt = CFL_MAX * min_dx / umax;
      if ( dt > max_dt)
         dt = max_dt;
      end
      t  = [ t; t(end) + dt ];
      fprintf( ['   New Time-Step: ', num2str(dt) '\n'] );

      if dt < 1.0e-10
         fprintf( 'Simulation has gone unstable. Exiting.\n' );
         return;
      end

   end

   %XXX debug.
   assignin( 'base', 'uxu', reshape( uxu, n * mz, n * mx, length(t) ) );
   assignin( 'base', 'uzu', reshape( uzu, n * mz, n * mx, length(t) ) );


end
