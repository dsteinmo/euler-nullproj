function [ N ] = compute_divergence_nullspace( D, n, mx, mz );
% [ N ] = compute_divergence_nullspace( D, n, mx, mz )
%
% Computes the nullspace basis of the discontinuous divergence
% operator element-by-element for speed.
%
%  Takes 4 arguments:
%
%     D         - divergence operator.
%     n, mx, mz - grid parameters.
%
%  Returns 1 argument:
%
%     N         - null space basis of D.
%
%  7 July 2015
%  Sumedh Joshi
%  Cornell University


   % Permute the divergence matrix into element-first indexing.
   r  = size( D, 1 );
   Dxe = D( :, 1:r );
   Dze = D( :, r+1:2*r );
   Dxe = permute_z2e( Dxe', n, mx, mz );
   Dxe = permute_z2e( Dxe', n, mx, mz );
   Dze = permute_z2e( Dze', n, mx, mz );
   Dze = permute_z2e( Dze', n, mx, mz );

   % Loop over elements, computing then null space basis within each.
   r = size( D, 2 );
   N = zeros( r, mx * mz * ( n^2 + 1) ); % Need to prove this is the dimension of N.
   %N = [];
   count = 1;
   for ii = 1:mx*mz

         ii
         iistart = n * n * ( ii - 1 ) + 1;
         iiend   = n * n * ( ii - 1 ) + n * n;

         % Compute the local null basis.
         iiD    = [ Dxe(iistart:iiend, iistart:iiend) Dze(iistart:iiend, iistart:iiend) ];
         iiNloc = null( full( iiD ) );

         % Split up the null spaces and put them into the right place in the global null matrix.
         iiN = zeros( r, size( iiNloc, 2 ) );
         iiN( iistart:iiend, : ) = iiNloc( 1:end/2, : );
         iistart = iistart + n * n * mx * mz;
         iiend   = iiend   + n * n * mx * mz;
         iiN( iistart:iiend, : ) = iiNloc( end/2 + 1:end, : );
         N(:, count:count+size(iiN,2) - 1 ) = iiN; % This seems wonky but its faster.
         count = count + size(iiN,2);
   end

   % Permute the null space basis into z-first indexing.
   r = n * n * mx * mz;
   Nxe = N( 1:r, : );
   Nze = N( r+1:end, : );
   Nxe = permute_e2z( Nxe', n, mx, mz );
   Nxe = permute_e2z( Nxe', n, mx, mz );
   Nze = permute_e2z( Nze', n, mx, mz );
   Nze = permute_e2z( Nze', n, mx, mz );
   N   = [Nxe; Nze];

end
