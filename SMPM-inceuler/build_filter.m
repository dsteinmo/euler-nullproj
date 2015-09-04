function exp_filter = build_filter( n, knots, weights, order )
% [exp_filter] = build_filter( n, knots, weights, order )
%
%  Builds an exponential filtering matrix on the quadrature weights/knots of order "order".
%
%  Takes 4 arguments:
%
%     n       - number of GLL knots.
%     knots   - GLL quadrature knots.
%     weights - GLL quadrature weights.
%     order   - exponential filter order.
%
%  Returns 1 argument:
%
%     exp_filter - an n x n matrix that is the filtering operator.
%
%
%  4 Sept 2015
%  Sumedh Joshi
%  Cornell University

  % Preallocate some variables.
  B = zeros( n, n );
  L = zeros( n, n );
  M = zeros( n, n );
  C = zeros( n, n );
  W = zeros( n, n );
  P = zeros( n, 1 );

  % B matrix -> Legendre polyn. evaluated in GLL knots
  for j = 1:n
     P(1) = 1.;
     P(2) = knots(j);
     for k = 2:n-1
        kp    = k + 1;
        km    = k - 1;
        P(kp) = (((2.0 * (k - 1)) + 1.) * P(2) * P(k) / (k + 1 - 1)) - ((k - 1) * P(km) / (k + 1 - 1));
     end

     for i = 1:n
        B(j, i) = P(i);
     end
  end

  for i = 1:n-1
     C(i, i) = (i - 1.) + 0.5;
     W(i, i) = weights(i);
  end

  C(n, n) = n / 2.;
  W(n, n) = weights(n);

  M = C * B';
  M = M * W;

  % Exponential filter
  alpha = -log10( 1.0e-16 );

  for k = 1:n
     L(k, k) = exp( -alpha * ((k - 1.0) / (n - 1.0)).^order );
  end

  exp_filter = B * L;
  exp_filter = exp_filter * M;

end
