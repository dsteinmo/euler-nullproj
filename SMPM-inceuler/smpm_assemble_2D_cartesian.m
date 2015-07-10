function [ D2x, D2z, E0, E1, B0, B1 ] = smpm_assemble_2D_cartesian( n, mx, mz, Lx, Lz )
% [Dx, Dz, E0, E1, B0, B1] = smpm_assemble_2D_cartesian( n, mx, mz, Lx, Lz )
%
%  Builds the domain decomposition of the spectral multidomain penalty method
%  discretization of the Poisson equation on a uniform cartesian grid.
%  Returns all matrices as sparse objects.
%
%  To assemble the SMPM matrix, compute:
%
%     SMPM = Dx * Dx + Dz * Dz + E;
%
%  Takes 5 arguments:
%
%    n  - Number of GLL points per direction, per subdomain.
%    mx - Number of subdomains in the x-direction.
%    mz - Number of subdomains in the z-direction.
%    Lx - Length between the first and last elements in the x-direction.
%    Lz - Length between the first and last elements in the z-direction.
%
%
%  24 June 2015
%  Sumedh Joshi

% Set the polynomial order in each element in each direction.
p = n - 1;

% Set the patching penalty factor.
F = 1.0;

% Compute the element-wise coordinates.

    % Compute the GLL gridpoints.
    n            = p;
    [x, w, junk] = lglnodes( p );
    x            = flipud( x )';

    % Compute the spectral differentiation matrix of order 2.
    D2 = zeros( n+1, n+1 );
    for ii = 1:n+1
        D2(:, ii) = make_lagrange( x, x, ii, 2 );
    end

    % Compute the spectral differentiation matrix of order 1.
    D1 = zeros( n+1, n+1 );
    for ii = 1:n+1
        D1(:, ii) = make_lagrange( x, x, ii, 1 );
    end
    x = x';

% Build the grid in xi-first indexing.
n      = n + 1;
r      = n * n * mx * mz;
%[x, z] = smpm_build_cartesian_mesh( n, mx, mz, [0, Lx], [0, Lz] );

% Set some deformation constants (in place of maps since we assume the grid is plaid).
hx = Lx / mx;
hz = Lz / mz;

% Build a 1D SMPM matrix in x.

   % The 1st derivative in x.
   Dx  = (2 / hx) * kron( speye( mx, mx ), D1 );

   % The 2nd derivative in x.
   Dxx = (2 / hx)^2 * kron( speye( mx, mx ), D2 );

   % The C0 continuity conditions.
   C0x = sparse( Dxx * 0 );
   for ii = 1:mx - 1
      left                        = ii * n;
      right                       = ii * n + 1;
      C0x(left:right, left:right) = [1 -1; -1 1];
   end

   % The C1 continuity conditions.
   C1x = sparse( Dxx * 0 );
   for ii = 1:mx - 1
      left                        = ii * n;
      right                       = ii * n + 1;
      C1x(left:right, left:right) = [ 1 -1; 1 -1];
   end

   % The Neumann boundary conditions.
   B1x           = sparse( Dxx * 0 );
   B1x(1,   1)   = -1;
   B1x(end, end) =  1;

   % Set the penalty parameter.
   w0     = w(1);
   taumin = (1 / w0) * (1 + 2*w0 - 2*sqrt( w0^2 + w0 )) * 2 / hx;
   taumax = (1 / w0) * (1 + 2*w0 + 2*sqrt( w0^2 + w0 )) * 2 / hx;
   tau    = (taumax + taumin) * 0.25;
   tau    = 1.0;

   % Assemble the matrix.
   Ex = tau * F * (C0x + sparse( C1x * Dx )) + tau * sparse( B1x * Dx );

% Build a 1D SMPM matrix in z.

   % The 1st derivative in z.
   Dz  = (2 / hz) * kron( speye( mz, mz ), D1 );

   % The 2nd derivative in z.
   Dzz = (2 / hz)^2 * kron( speye( mz, mz ), D2 );

   % The C0 continuity conditions.
   C0z = sparse( Dzz * 0 );
   for ii = 1:mz - 1
      left                        = ii * n;
      right                       = ii * n + 1;
      C0z(left:right, left:right) = [1 -1; -1 1];
   end

   % The C1 continuity conditions.
   C1z = sparse( Dzz * 0 );
   for ii = 1:mz - 1
      left                        = ii * n;
      right                       = ii * n + 1;
      C1z(left:right, left:right) = [ 1 -1; 1 -1];
   end

   % The Neumann boundary conditions.
   B1z           = sparse( Dzz * 0 );
   B1z(1,   1)   = -1;
   B1z(end, end) =  1;

   % Set the penalty parameter.
   w0     = w(1);
   taumin = (1 / w0) * (1 + 2*w0 - 2*sqrt( w0^2 + w0 )) * 2 / hz;
   taumax = (1 / w0) * (1 + 2*w0 + 2*sqrt( w0^2 + w0 )) * 2 / hz;
   tau    = (taumax + taumin) * 0.25;
   tau    = 1.0;

   % Assemble the matrix.
   Ez = tau * F * (C0z + sparse( C1z * Dz )) + tau * sparse( B1z * Dz );

% Build the 2D non-divergent part of the operator.
E    = kron( speye( mx * n ), ...
             sparse( Ez ) ) + ...
       kron( sparse( Ex ), ...
             speye( mz * n ) );

E0   = kron( speye( mx * n ), sparse( C0z ) ) + kron( sparse( C0x ), speye( mz * n ) );
E1   = kron( speye( mx * n ), sparse( C1z * Dz ) ) + kron( sparse( C1x * Dx ), speye( mz * n ) );
B1   = kron( speye( mx * n ), sparse( B1z * Dz ) ) + kron( sparse( B1x * Dx ), speye( mz * n ) );
B0   = kron( speye( mx * n ), sparse( abs(B1z) ) ) + kron( sparse( abs(B1x) ), speye( mz * n ) );

% Build two-dimensional versions of the derivative operator.
D2z = kron( eye( mx * n ), Dz );
D2x = kron( Dx, eye( mx * n ) );

end
