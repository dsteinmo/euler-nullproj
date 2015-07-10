%
% Computes permutation of grid functions from element-based to z-first indexing.
function A = permute_e2z( A, n, mx, mz )

   counter = 1;
   swap = 0 * A;
   for xelt = 1:mx
      for jj = 1:n
         for zelt = 1:mz
            for ii = 1:n

               ndx_elt = ( xelt - 1 ) * n * n + ( zelt - 1 ) * n * n * mx + ( ii - 1 ) * n + jj;
               ndx_z   = counter;
               counter = counter + 1;
               swap( ndx_z, : ) = A( ndx_elt, : );

            end
         end
      end
   end
   A = swap;
