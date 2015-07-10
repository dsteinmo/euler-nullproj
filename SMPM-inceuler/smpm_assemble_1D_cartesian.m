function [D, G, E] = smpm_assemble_1D_cartesian( n, mx, Lx )
% [D, G, E] = smpm_assemble_1D_cartesian( n, mx, Lx )
%
%  To assemble the SMPM matrix, compute:
%
%     SMPM = D*G + E;
%
%  Takes 3 arguments:
%
%    n  - Number of GLL points per direction, per subdomain.
%    mx - Number of subdomains in the x-direction.
%    Lx - Length between the first and last elements in the x-direction.
%
%  Returns 4 values:
%
%    A -
%    E -
%    B -
%    S -
%
%  24 June 2015.
%  Sumedh Joshi

% Parse the arguments.

% Set the polynomial order in each element in each direction.
p = n - 1;

% Set the patching penalty factor.
F = 1000;

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

% Re-set the value of n.
n = n + 1;

% Set some deformation constants (in place of maps since we assume the grid is plaid).
hx = Lx / mx;

% Build the D and G matrices (in 1D they are the same).
D = ( 2 / hx ) * kron( eye( mx, mx ), D1 );
G = D;
E = 0 * G;

% Build a 1D SMPM matrix in x.

   % The 1st derivative in x.
   Dx  = (2 / hx) * kron( speye( mx, mx ), D1 );

   % The 2nd derivative in x.
   Dxx = (2 / hx)^2 * kron( speye( mx, mx ), D2 );

   % The C0 continuity conditions.
   C0x = 0 * E;
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

   % Assemble the matrix.
   %SMPMx = Dxx + tau * F * (C0x + sparse( C1x * Dx )) + tau * sparse( B1x * Dx );

   % Assemble the E matrix.
   E = tau * F * (C0x + sparse( C1x * Dx )) + tau * sparse( B1x * Dx );

end
