function [ux uz] = apply_smpm_boundary_conditions( ux, uz, ux_b, uz_b, n, mx, mz, Lx, Lz, nu )
% [ux uz] = apply_smpm_boundary_conditions( ux, uz, ux_b, uz_b, n, mx, mz, Lx, Lz, boundaries, nu )
%
%  Sets the SMPM viscous boundary values for the RHS of the viscous equations.
%  Only works for Dirichlet conditions.
%
%  Takes 8 arguments:
%
%     ux, uz     - velocity vectors.
%     ux_b, uz_b - boundary values.
%     n, mx, mz  - grid discretization constants.
%     Lx, Lz     - penalty coefficients in the x- and z-direction.
%     boundaries - BC vector.
%     nu         - viscosity.
%
%  Returns 2 arguments:
%
%     ux, uz      - velocity.
%
%  27 August 2015
%  Sumedh Joshi
%  Cornell University

   % Set some constants.
   pd = n - 1;
   omega = 2.0 / ( pd * ( pd + 1.0 ) );
   tau_d = nu / 1.0 / omega^2;
   tau_n = nu / omega;
   rpk   = n * n * mx * mz;
   hx    = Lx / mx;
   hz    = Lz / mz;

   % XXX: only works for non-homogenous Dirichlet boundary conditions.

   % Apply the boundary condition to the left boundary.
   a      = 1
   z      = n * mx
   ux(a:z) = ux(a:z) + tau_d * ( 2 / hx ) * ux_b(a:z)

   % Apply the boundary condition to the right boundary.
   a      = rpk - (n * mx - 1)
   z      = rpk - 0
   ux(a:z) = ux(a:z) + tau_d * ( 2 / hx ) * ux_b(a:z)

   % Apply the boundary condition to the top boundary.
   a          = n * mx
   z          = rpk
   gap        = n * mx
   uz(a:gap:z) = uz(a:gap:z) + tau_d * ( 2 / hz ) * uz_b(a:gap:z)

   % Apply the boundary condition to the bottom boundary.
   a          = 1
   z          = rpk
   gap        = n * mx
   uz(a:gap:z) = uz(a:gap:z) + tau_d * ( 2 / hz ) * uz_b(a:gap:z)

end
