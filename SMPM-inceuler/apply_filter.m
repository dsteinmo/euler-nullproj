function Fq = apply_filter( n, mx, mz, q, filter )
% Fq = apply_filter( q, filter )
%
%  Applies an exponential filter as built by build_filter.
%
%  Takes 2 arguments:
%
%     q      - function to be filtered.
%     filter - filtering matrix.
%
%  Returns 1 argument:
%
%     Fq - filtered function.
%
%
%  4 Sept 2015
%  Sumedh Joshi
%  Cornell University

  % Permute into element-first indexing.
  q = permute_z2e( q, n, mx, mz );

  rhoN = zeros( n * n * mx * mz, 1 );
  for k = 0:mx*mz - 1
     rhoT  = zeros( n, n );
     rhoF1 = zeros( n, n );;

     % Extracting the subdomain information from global solution
     for i = 0:n-1
        for j = 1:n
           c            = (n * n * k) + (i * n) + j;
           rhoT(i+1, j) = q(c);
        end
     end

     % Now, filtering the solution at the subdomain
     rhoF1 = filter * rhoT;
     rhoF2 = rhoF1 * filter';

     % Returning the filtered solution to the global system
     for i = 0:n - 1
        for j = 1:n
           c       = (n * n * k) + (i * n) + j;
           rhoN(c) = rhoN(c) + rhoF2(i + 1, j);
        end
     end
  end

  q = rhoN;

  % Permute into xi-first indexing.
  Fq = permute_e2z( q, n, mx, mz );

end
