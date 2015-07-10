function [ vx, vz ] = apply_nullspace_projection( ux, uz, N, U, S, V );
% [vx, vz] = apply_nullspace_projection( ux, uz, N, U, S, V )
%
%  Assembles and solves a regularized least squares problem with a truncated
%  SVD.  The regularization parameter, tau, must be zero or positive for there
%  to exist a solution.  The resulting velocity vectors vx and vz are
%  exactly divergence-free.
%
%  Takes 5 arguments:
%
%     ux, uz  - velocity vectors.
%     N       - basis for the null space of the divergence operator.
%     U, S, V - truncated SVD as returned by setup_nullspace_projection().
%
%  Returns 2 arguments:
%
%     vx, vz - projection of (ux,uz) onto the null space of the divergence operator.
%
%  5 July 2015
%  Sumedh Joshi
%  Cornell University

   % Solve the problem with a truncated SVD.
   lambda = V * ( S \ ( U' * ( N' * [ux; uz] ) ) );
   v = N * lambda;
   vx = v(1:end/2);
   vz = v(end/2+1:end);

   % Check the accuracy with which we solved the problem.
   err = norm( U * ( S * ( V' * lambda ) )- N'*[ux; uz] );
   fprintf( ['   Error in projection solve:', num2str( err ) '\n'] );

end
