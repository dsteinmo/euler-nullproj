%
% Script to post-process/visualize solutions of the three methods for solving the
% Rayleigh-Taylor problem.

   % Name the three methods.
   method = { 'Poisson', 'Post-Projection', 'Penalized Projection', 'Post-Weak-Projection' };

   % Set some filenames in order.
   fnames = { 'djl_poisson.mat', 'djl_postproject.mat', 'djl_postnormal.mat', 'djl_postnull.mat' };
   fnames = { 'djl_poisson_hires.mat', 'djl_postproject_hires.mat', 'djl_postnormal_hires.mat', 'djl_postnull_hires.mat' };

   % Set the visualization window in x.
   win = [ 5, 10 ];

   % Set the time you want to extract.
   tplot = 48.0;

   % Figure 0: total number of time-steps.
   figure;
   hold all;
   for iiviz = 1:length(method)

      load( fnames{iiviz} );
      plot(t);

   end
   xlabel( 'time-step' );
   ylabel( 'physical time' );
   box on;
   legend( method, 'Location', 'NorthEast' );
   print_graphics( gcf, 'compare_djl_time', 1, 0, 0, 0 );


   % Figure 1: quiver plot of the final velocity.
   figure;
   for iiviz = 1:length(method)

      load( fnames{iiviz} );
      x = reshape( x, n * mz, n * mx );
      z = reshape( z, n * mz, n * mx );
      [junk ndx] = min( abs( t - tplot ) );
      subplot( 2, 2, iiviz );
      quiver( x, z, ux(:,:,ndx), uz(:,:,ndx), 1.0 );
      set( gca, 'xlim', win );
      title( [ method{iiviz} ': velocity' ] );

   end
   set( gcf,'pos',[2     5   948   690]);
   print_graphics( gcf, 'compare_djl_velocity', 1, 0, 0, 0 );

   % Figure 2: mass conservation.
   figure;
   hold all;
   for iiviz = 1:length(method)

      load( fnames{iiviz} );
      [ junk ndx ] = min( abs( t - tplot ) );
      mass = integrate_grid_function( reshape( rho, n * n * mx * mz, length(t) ), x, z, n, mx, mz );
      plot( t(1:ndx), mass(1:ndx) );

   end
   xlabel( 'time' );
   ylabel( 'total mass' );
   box on;
   legend( method, 'Location', 'NorthEast' );
   print_graphics( gcf, 'compare_djl_mass', 1, 0, 0, 0 );

   % Figure 3: continuity.
   figure;
   for iiviz = 1:length(method)

      load( fnames{iiviz} );
      ux = reshape( ux, r, length(t) );
      uz = reshape( uz, r, length(t) );
      [junk ndx] = min( abs( t - tplot ) );
      c = sqrt( sum( ( E0 * ux(:,1:ndx) ).^2 + ( E0 * uz(:,1:ndx) ).^2 ) );
      semilogy( t(1:ndx), c(1:ndx) );
      hold all;

   end
   xlabel( 'time' );
   ylabel( 'C^0 discontinuity norm' );
   legend( method , 'Location', 'NorthEast' );
   box on;
   print_graphics( gcf, 'compare_djl_continuity', 1, 0, 0, 0 );

   % Figure 4: divergence.
   figure;
   for iiviz = 1:length(method)

      load( fnames{iiviz} );
      [ junk ndx ] = min( abs( t - tplot ) );
      ux = reshape( ux, r, length(t) );
      uz = reshape( uz, r, length(t) );
      divu = sqrt( sum( ( D * [ux(:,1:ndx); uz(:,1:ndx)]).^2 ) );
      semilogy( t(1:ndx), divu );
      hold all;

   end
   xlabel( 'time' );
   ylabel( 'norm of divergence' );
   legend( method , 'Location', 'NorthEast' );
   box on;
   print_graphics( gcf, 'compare_djl_divergence', 1, 0, 0, 0 );

   % Figure 5: contour plot of the final density
   figure;
   for iiviz = 1:length(method)

      [ junk, ndx ] = min( abs( t - tplot ) );

      load( fnames{iiviz} );
      subplot( 2, 2, iiviz );
      [ junk ndx ] = min( abs( t - tplot ) );
      rhob = reshape( rhob, n * mz , n * mx );
      x = reshape( x, n * mz, n * mx );
      z = reshape( z, n * mz, n * mx );
      contour( x, z, rhob + rho(:,:,ndx), 5 );
      set( gca, 'xlim', win );
      title( [ method{iiviz} ': density' ] );
      colorbar;

   end
   set( gcf,'pos',[2     5   948   690]);
   print_graphics( gcf, 'compare_djl_density', 1, 0, 0, 0 );

   % Figure 5b: contour plot of the final u-velocity.
   figure;
   for iiviz = 1:length(method)

      [ junk, ndx ] = min( abs( t - tplot ) );

      load( fnames{iiviz} );
      subplot( 2, 2, iiviz );
      [ junk ndx ] = min( abs( t - tplot ) );
      rhob = reshape( rhob, n * mz , n * mx );
      x = reshape( x, n * mz, n * mx );
      z = reshape( z, n * mz, n * mx );
      ux = reshape( ux, n * mz, n * mx, length(t) );
      contour( x, z, ux(:,:,ndx), 5 );
      set( gca, 'xlim', win );
      title( [ method{iiviz} ': u-velocity' ] );
      colorbar;

   end
   set( gcf,'pos',[2     5   948   690]);
   print_graphics( gcf, 'compare_djl_ux', 1, 0, 0, 0 );

   % Figure 5c: contour plot of the final w-velocity.
   figure;
   for iiviz = 1:length(method)

      [ junk, ndx ] = min( abs( t - tplot ) );

      load( fnames{iiviz} );
      subplot( 2, 2, iiviz );
      [ junk ndx ] = min( abs( t - tplot ) );
      rhob = reshape( rhob, n * mz , n * mx );
      x = reshape( x, n * mz, n * mx );
      z = reshape( z, n * mz, n * mx );
      uz = reshape( uz, n * mz, n * mx, length(t) );
      contour( x, z, uz(:,:,ndx), 5 );
      set( gca, 'xlim', win );
      title( [ method{iiviz} ': w-velocity' ] );
      colorbar;

   end
   set( gcf,'pos',[2     5   948   690]);
   print_graphics( gcf, 'compare_djl_uz', 1, 0, 0, 0 );

   % Figure 6: show norm of velocity field to get a sense of instability.
   figure;
   for iiviz = 1:length(method)

      load( fnames{iiviz} );
      ux = reshape( ux, r, length(t) );
      uz = reshape( uz, r, length(t) );
      normu = sqrt( sum( ux.^2 + uz.^2 ) );
      semilogy( t, normu );
      hold all;

   end
   xlabel( 'time' );
   ylabel( 'norm of velocity' );
   legend( method , 'Location', 'NorthEast' );
   box on;
   print_graphics( gcf, 'compare_djl_velocity_norm', 1, 0, 0, 0 );

   % Figure 7: kinetic energy conservation.
   figure;
   hold all;
   for iiviz = 1:length(method)

      load( fnames{iiviz} );
      [ junk ndx ] = min( abs( t - tplot ) );
      u = sqrt( ux.^2 + uz.^2 );
      KE = 0 * u;
      rhob = reshape( rhob, n * mz, n * mx );
      for ii = 1:length(t)
         KE(:,:,ii) = 0.5 * ( rho0 + rhob + rho(:,:,ii) ).* u(:,:,ii);
      end
      energy = integrate_grid_function( reshape( KE, n * n * mx * mz, length(t) ), x, z, n, mx, mz );
      plot( t(1:ndx), energy(1:ndx) );

   end
   xlabel( 'time' );
   ylabel( 'total kinetic energy' );
   box on;
   legend( method, 'Location', 'NorthEast' );
   print_graphics( gcf, 'compare_djl_energy', 1, 0, 0, 0 );
