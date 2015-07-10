function I = integrate_grid_function( u, x, z, n, mx, mz );
% I = integrate_grid_function( u, x, z, n, mx, mz );
%
%  Computes the GLL quadrature integral of the grid function u
%  defined on a cartesian multi-element grid.
%
%  Takes 7 arguments:
%
%     u          - grid function to be integrated.
%     x, z       - grid coordinates.
%     n, mx, mz  - grid discretization constants.
%
%  Returns 1 arguments:
%
%     I          - value of the integral.
%
%  7 July 2015
%  Sumedh Joshi
%  Cornell University

   % Get some quadrature weights.
   [junk w P] = lglnodes( n - 1 );

   % Get some element widths.
   Lx = max(x) - min(x);
   Lz = max(z) - min(z);
   hx = Lx / mx;
   hz = Lz / mz;

   % Loop through grid function, mutiplying by quadrature weights.
   for ii = 1:mx
      for jj = 1:n
         for kk = 1:mz
            for ll = 1:n

               ndx = ( ii - 1 ) * n * n * mz + ( jj  - 1 ) * n * mz +  ( kk - 1 ) * n + ll;
               u(ndx,:) = w(ll) * w(jj) * u(ndx,:) * ( hx / 2 ) * ( hz / 2 );

            end
         end
      end
   end

   % Sum the kernel.
   I = sum(u);

end
