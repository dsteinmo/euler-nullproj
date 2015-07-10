function Aq = apply_smpm_advection( q, Dx, Dz, ux, uz, Lx, Lz, n, mx, mz );
% Aq = apply_smpm_advection( q, Dx, Dz, ux, uz, Lx, Lz, n, mx, mz )
%
%  Applies the SMPM advection operator on a cartesian, Lx x Lz domain with a
%  cartesian mx x mz grid with n GLL points per element.
%
%  Takes 7 arguments:
%
%     q          - field variable to be advected.
%     Dx, Dz     - multi-element multi-dimensional spectral differentiation matrices.
%     ux, uz     - velocity vectors.
%     Lx, Lz     - penalty coefficients in the x- and z-direction.
%     n, mx, mz  - grid discretization constants.
%
%  Returns 1 arguments:
%
%     Aq           - advected quantity.
%
%  7 July 2015
%  Sumedh Joshi
%  Cornell University

   % Set some constants for penalty parameters.
   pdegree = n - 1;
   omega   = 2 / ( pdegree * ( pdegree + 1 ) );
   hx      = Lx / mx;
   hz      = Lz / mz;
   tau     = 1 / omega;

   % Compute the gradient of the field.
   q_x = Dx * q;
   q_z = Dz * q;

   % Compute the divergence of the velocity times the field.
   div_uq = Dx * ( ux .* q ) + Dz * ( uz .* q );

   % Compute the advection operator in skew-symmetric form.
   Aq = -1.0 * ( 0.5 * ux .* q_x + uz .* q_z ) - 0.5 * div_uq;

   % Compute u-dot-n in the vertical direction (normal vectors always point up).
   udotn = uz;

   % Apply inter-element patching conditions along horizontal interfaces.
   for ii = 1: n * mx

      % Loop over the horizontal interfaces.
      for jj = 1: mz - 1

         % Get the indices for the top/bottom of this interface.
         bottom = ( ii - 1 ) * ( mz * n ) + jj * n;
         top    = ( ii - 1 ) * ( mz * n ) + jj * n + 1;

         % Fix the fact that my normal vectors always point up.
         udotn(top) = -1 * udotn(top);

         % Check for inflow below the interface.
         if ( udotn(bottom) < 0. )

            % Apply flux to the bottom element.
            Aq(bottom) = Aq(bottom) - tau * ( 2 / hz ) * abs( udotn(bottom) ) * ( q(bottom) - q(top) );

         else

            % Apply flux to the top element.
            Aq(top) = Aq(top) - tau * ( 2 / hz ) * abs( udotn(top) ) * ( q(top) - q(bottom) );

         end
      end
   end

   % Compute u-dot-n in the horizontal (normal vectors always point left).
   udotn = -ux;

   % Apply inter-element patching conditions along vertical interfaces.
   for ii = 1: mx - 1

      % Loop over the nodes in this vertical interface.
      for jj = 1: n * mz

         % Get the indices for the left and right sides of this interface.
         left  = ii * n * n * mz - n * mz + jj;
         right = ii * n * n * mz + jj;

         % Correct the normal vector to point outward on the left.
         udotn(left) = -1 * udotn(left);

         % Check for inflow on the left side of the interface.
         if ( udotn(left) < 0.0 )

            % Apply flux to the left.
            Aq(left) = Aq(left) - tau * ( 2 / hx ) * abs( udotn(left) ) * ( q(left) - q(right) );

         else

            % Apply flux to the right.
            Aq(right) = Aq(right) - tau * ( 2 / hx ) * abs( udotn(right) ) * ( q(right) - q(left) );

         end
      end
   end

end
