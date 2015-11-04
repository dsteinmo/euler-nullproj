%%
%
% A matrix-storage implementation of GMRES with Householder reflections via
% the algorithm proposed in Walker 1988.  A is either a square matrix or a
% struct that has .handle as a function handle and can be called by
% inputs.handle( inputs ) and that returns y = A(x), where y = A*x.
%
% [x err m] = compute_gmres_householder(A,b,x0,TOL,MAXIT);
%
% - OR -
%
% [x err m] = compute_gmres_householder(A,b,x0,TOL,MAXIT,PRE);
%
% PRE is either a (right) preconditioning matrix, or a cell array of
% additive right preconditioning matrices.
%
% err is the history of the error over m iterations, and x is the solution.
%
% 27 Oct 2013
% Sumedh Joshi
% Cornell University

function [x err_array mout] = compute_gmres_householder(varargin)
% function [x err mout] = compute_gmres_householder(A,b,x0,TOL,MAXIT,PRE)

    %
    % Argument parsing.
    if nargin < 6
        A     = varargin{1};
        b     = varargin{2};
        x0    = varargin{3};
        TOL   = varargin{4};
        MAXIT = varargin{5};
        PRE   = speye(size(A));
    else
        A     = varargin{1};
        b     = varargin{2};
        x0    = varargin{3};
        TOL   = varargin{4};
        MAXIT = varargin{5};
        PRE   = varargin{6};
    end

    %
    % See if the first input was a struct.
    if isstruct( A )
       inputs = A;
       A = inputs.handle;
    end

    %
    % Get some constants.
    n = size(x0,1);

    %
    % Set the maxit.
    MAXIT = min([n,MAXIT]);

    %
    % Adjust the tolerance to be a relative error only if the norm of the right-hand-side is bigger than one.
    % This way, we guarantee that the relative error is always less than |b|.
    NORM_RHS = max( 1.0, norm(b) );

    %
    % Going forward, the comments and variable names
    % correspond exactly to the steps outlined
    % in Walker's paper.
    %
    %   Except for G, my storage for the Givens rotation matrices; the
    %   paper uses the letter J to denote the Givens rotations.

    %
    % 0. The outer iteration (i.e. wrapping the restarts) begin here.
    converged = 0;
    while ~converged

        %
        % Set up some variables I'll need.
        P = [];
        G = [];
            %
            % G is the givens rotation storage vector.  The Givens rotation
            % corresponding to G(:,j) is
            %
            % | G(1,j)  G(2,j) |
            % |-G(2,j)  G(1,j) |

        %
        % 1. Compute r0 and find the appropriate Householder reflector.
        if isa(A,'function_handle')
            if iscell( PRE )
               inputs.x = 0;
               for iipre = 1:length(PRE)
                  inputs.x = inputs.x + PRE{iipre} \ x0;
               end
            else
               inputs.x = PRE\x0;
            end
            r0 = b - A(inputs);
        else
            if iscell( PRE )
               applypre = 0;
               for iipre = 1:length(PRE)
                  applypre = applypre + PRE{iipre} \ x0;
               end
            else
               applypre = PRE\x0;
            end
            r0 = b - A*(applypre);
        end
        P(:,1) = get_householder(r0,1);
        w      = apply_householder(r0,P(:,1));

        %
        % Set the initial error in the array of errors.
        err_array    = [];
        err_array(1) = norm( r0 );
        err_array(1) = w(1);

        %
        % 2. For m=1,2,...MAXIT do:
        m = 1;
        while(1)

            %
            % a) Evalute v = PmPm-1...P1AP1...Pmem.
            v = zeros(n,1);
            v(m) = 1;
            for jj = m:-1:1
                v = apply_householder(v,P(:,jj));
            end
            applypre = 0.0;
            if iscell(PRE)
               for iipre = 1:length(PRE)
                  applypre = applypre + PRE{iipre} \ v;
               end
            else
               applypre = PRE \ v;
            if isa(A,'function_handle')
                inputs.x = applypre;
                v = A(inputs);
            else
                v = A*(applypre);
            end
            for jj = 1:m
                v = apply_householder(v,P(:,jj));
            end

            %
            % b) if v(m+1) = ... = v(n) = 0, then proceed to step (e),
            % otherwise continue.
            if ~all(v(m+1:end) <= eps)

                %
                % c) Determine Pm+1 with a Householder vector having first m
                % components 0, such that Pm+1v has zero components after the
                % (m+1)st.
                smallP = get_householder(v(m+1:end),1);
                P(:,m+1) = [zeros(m,1); smallP];

                %
                % d) Overwrite v <-- Pm+1v. 
                v = apply_householder(v,P(:,m+1)); 
            else
                P(:,m+1) = zeros(n,1);
            end

            %
            % e) If m > 1, overwrite v <-- Jm-1..J1v. 
            if m > 1
                for jj = 1:m-1

                    %
                    % Construct this Givens rotation. 
                    J(1,1) =  G(1,jj);
                    J(1,2) =  G(2,jj); 
                    J(2,1) = -G(2,jj);
                    J(2,2) =  G(1,jj); 

                    %
                    % Apply this Givens rotation. 
                    v(jj:jj+1) = J*v(jj:jj+1); 

                end
            end

            %
            % f) If v(m+1) = 0, proceed to step (i), otherwise, continue. 
            if v(m+1) ~= 0 

                %
                % g) Determine Jm acting on components m and m+1 such that 
                % (Jmv)(m+1) is zero. 
                thisgivens = planerot([v(m),v(m+1)]'); 
                G(:,m) = [thisgivens(1,1), thisgivens(1,2)]'; 

                %
                % h) Overwrite v <-- Jmv and w <-- Jmw.
                v(m:m+1) = thisgivens*v(m:m+1);
                w(m:m+1) = thisgivens*w(m:m+1);

            end

            %
            % i) Set the matrix Rm = ... etc.
            if m == 1
                R = v;
            else
                R = [R v];
            end

            %
            % j) If |w(m+1)| < TOL or m = MAXIT, then solve for ym and
            % overwrite x0 with xm; otherwise increment m.
            err_array( m + 1 ) = w(m+1);
            if abs(w(m+1)) < TOL * NORM_RHS || m == MAXIT

                %
                % 3. Solve for ym and overwrite x0 <--- xm.

                    %
                    % a) Determine ym which minimizes the least-squares problem
                    % for ||w - Rm y|| by solving the m x m upper triangular
                    % system with the first m rows of Rm as the coefficient
                    % matrix and the first m components of w as the right hand
                    % side.
                    Rm = R(1:m,1:m);
                    wm  = w(1:m);

                    %
                    % Compute the least squares solution via
                    % backsubstitution.
                    sm    = zeros(size(wm));
                    sm(m) = wm(m)/Rm(m,m);
                    for ii = m-1:-1:1
                       sm(ii) = (wm(ii) - (Rm(ii,ii+1:m)*sm(ii+1:m)))/Rm(ii,ii);
                    end
                    ym = sm;

                    %
                    % b) For k = 1,...,m do:
                    %   overwrite x0 <-- x0 + ym(k)P1...Pkek.
                    for k = 1:m
                        ek      = zeros(n,1);
                        ek(k)   = 1;
                        xk      = ek;
                        for jj = k:-1:1
                            xk = apply_householder(xk,P(:,jj));
                        end
                        x0 = x0 + ym(k)*xk;
                    end

                    %
                    % c) If |w(m+1)| < TOL accept x0 as the solution; otherwise
                    % set m = 1 and start over.
                    if abs(w(m+1)) < TOL * NORM_RHS
                        converged = 1;
                        if iscell( PRE )
                           x = 0;
                           for iipre = 1:length(PRE)
                              x = x + PRE{iipre} \ x0;
                           end
                        else
                           x = PRE\x0;
                        end
                        mout = m;
                        err = w(m+1);
                        err = err_array;
                        break;
                    else
                        mout = m;
                        err = w(m+1);
                        err = err_array;
                        m = 1;
                        break;
                    end
            else

                %
                % Increment the counter.
                m = m + 1;

            end

        end

    end

end
