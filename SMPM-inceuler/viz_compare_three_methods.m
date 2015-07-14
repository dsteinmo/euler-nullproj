%
% Script to post-process/visualize solutions of the three methods for solving the
% Rayleigh-Taylor problem.

   % Name the three methods.
   method = { 'Poisson', 'Post-Projection', 'Nullspace-Projection', 'Post-Weak-Projection' };

   % Set some filenames in order.
   fnames = { 'raytay_poisson.mat', 'raytay_postproject.mat', 'raytay_nullspace.mat', 'raytay_postnull.mat' };

   % Set the visualization window.
   win = [ 1 / 2 - 0.2 * 1,  1 / 2 + 0.2 * 1 ];
%   win = [ 0, 1 ];

   % Set the time you want to extract.
   tplot = 0.85;

   % Figure 1: quiver plot of the final velocity.
   figure;
   for iiviz = 1:length(method)


      load( fnames{iiviz} );
      [junk ndx] = min( abs( t - tplot ) );
      subplot( 2, 2, iiviz );
      quiver( x, z, ux(:,ndx), uz(:,ndx), 5.0 );
      set( gca, 'xlim', win, 'ylim', win );
      title( [ method{iiviz} ': velocity' ] );

   end
   set( gcf,'pos',[2     5   948   690]);
   print_graphics( gcf, 'compare_velocity', 1, 0, 0, 0 );

   % Figure 2: mass conservation.
   figure;
   hold all;
   for iiviz = 1:length(method)

      load( fnames{iiviz} );
      [ junk ndx ] = min( abs( t - tplot ) );
      mass = integrate_grid_function( rho, x, z, n, mx, mz );
      plot( t(1:ndx), mass(1:ndx) );

   end
   xlabel( 'time' );
   ylabel( 'total mass' );
   box on;
   legend( method, 'Location', 'SouthEast' );
   print_graphics( gcf, 'compare_mass', 1, 0, 0, 0 );

   % Figure 3: continuity.
   figure;
   for iiviz = 1:length(method)

      load( fnames{iiviz} );
      [junk ndx] = min( abs( t - tplot ) );
      c = sqrt( sum( ( E0 * ux(:,1:ndx) ).^2 + ( E0 * uz(:,1:ndx) ).^2 ) );
      semilogy( t(1:ndx), c(1:ndx) );
      hold all;

   end
   xlabel( 'time' );
   ylabel( 'C^0 discontinuity norm' );
   legend( method , 'Location', 'SouthEast' );
   box on;
   print_graphics( gcf, 'compare_continuity', 1, 0, 0, 0 );

   % Figure 4: divergence.
   figure;
   for iiviz = 1:length(method)

      load( fnames{iiviz} );
      [ junk ndx ] = min( abs( t - tplot ) );
      divu = sqrt( sum( ( D * [ux(:,1:ndx); uz(:,1:ndx)]).^2 ) );
      semilogy( t(1:ndx), divu );
      hold all;

   end
   xlabel( 'time' );
   ylabel( 'norm of divergence' );
   legend( method , 'Location', 'SouthEast' );
   box on;
   print_graphics( gcf, 'compare_divergence', 1, 0, 0, 0 );

   % Figure 5: contour plot of the final density
   figure;
   for iiviz = 1:length(method)

      [ junk, ndx ] = min( abs( t - tplot ) );

      load( fnames{iiviz} );
      subplot( 2, 2, iiviz );
      [ junk ndx ] = min( abs( t - tplot ) );
      contourf( reshape(x, n * mx, n * mz ), reshape( z, n * mx, n * mz  ), 1000.0 + reshape( rho(:,ndx), n * mx, n * mz ), 25 )
      set( gca, 'xlim', win, 'ylim', win );
      title( [ method{iiviz} ': density' ] );
      colorbar;

   end
   set( gcf,'pos',[2     5   948   690]);
   print_graphics( gcf, 'compare_density', 1, 0, 0, 0 );

   % Figure 6: show norm of velocity field to get a sense of instability.
   figure;
   for iiviz = 1:length(method)

      load( fnames{iiviz} );
      normu = sqrt( sum( ux.^2 + uz.^2 ) );
      semilogy( t, normu );
      hold all;

   end
   xlabel( 'time' );
   ylabel( 'norm of velocity' );
   legend( method , 'Location', 'SouthEast' );
   box on;
   print_graphics( gcf, 'compare_velocity_norm', 1, 0, 0, 0 );
