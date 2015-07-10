function [ U S V ] = setup_nullspace_projection( tau, N, E );
% [ U S V ] = setup_nullspace_projection( tau, N, E )
%
% Computes and truncates the SVD of the normal equations given by
% I + tau * ( E * N )' * (E * N ).
% The outputs of this routine are used in solving the normal
% equations.
%
%  Takes 3 arguments:
%
%     tau    - regularization coefficient, positive.
%     N      - null space basis of the divergence operator.
%     E      - inter-element continuity and boundary condition operator.
%
%  Returns 3 arguments:
%
%     U, S, V - Truncated SVD of the normal equations.
%
%  7 July 2015
%  Sumedh Joshi
%  Cornell University

   % Set up the normal equations.
   r = size( N, 2 );
   %T = eye( r, r ) + tau * ( E * N )' * E * N;
   T = N'*N + tau * ( E * N )' * E * N;

   % Solve the problem with a truncated SVD.
   [U S V] = svd( T );
   rT = rank( T );
   U = U( :   , 1:rT);
   S = S( 1:rT, 1:rT);
   V = V( :   , 1:rT );

end
