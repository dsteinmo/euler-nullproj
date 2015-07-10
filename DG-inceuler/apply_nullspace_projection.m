function [ vx, vz ] = apply_nullspace_projection( ux, uz, tau, N, E );
% [vx, vz] = apply_nullspace_projection( )
%
%  Assembles and solves a regularized least squares problem with a truncated
%  SVD.  The regularization parameter, tau, must be zero or positive for there
%  to exist a solution.  The resulting velocity vectors vx and vz are
%  exactly divergence-free.
%
%  Takes 5 arguments:
%
%     ux, uz - velocity vectors.
%     tau    - regularization coefficient, positive.
%     N      - null space basis of the divergence operator.
%     E      - inter-element continuity and boundary condition operator.
%
%  Returns 2 arguments:
%
%     vx, vz - projection of (ux,uz) onto the null space of the divergence operator.
%
%  5 July 2015
%  Sumedh Joshi
%  Cornell University

   % Set up the normal equations.
   r = size( N, 2 );
   T = eye( r, r ) + tau * ( E * N )' * E * N;

   % Solve the problem with a truncated SVD.
   [U S V] = svd( T );
   rT = rank( T );
   lambda =  V(:,1:rT) * ( S(1:rT,1:rT) \ ( U(:,1:rT)' * ( N' * [ ux; uz ] ) ) );
   v = N * lambda;
   vx = v(1:end/2);
   vz = v(end/2+1:end);

end
