%consider ODE y' = f(y(t),t)
%as input, takes:
     %tn(1) = tn+1 (next time)
     %tn(2) = tn   (current time)
     %tn(3) = tn-1 (previous time)
     %tn(4) = tn-2 (two times ago)
     %bn(1) = coefficient of f(tn,yn)
     %bn(2) = coefficient of f(tn-1,yn-1)
     %bn(3) = coefficient of f(tn-2,yn-2) 
     %tn(2) to tn(1) are integration limits
     %the right-hand side is interpolated from tn to tn+1
     %using tn,tn-1, and tn-2
     %NOTE: ode is evolved afterwards via
     %yn+1 = yn + bn(1)*f(tn,yn)+bn(2)*f(tn-1,yn-1)+bn(3)*f(tn-2,yn-2)
     
     %standard results: [3 2 1 0] -> gives AB3 coefficients without dt factor
     %                  [2 1 0 0] -> gives AB2 coefficients without dt factor
     %                  [1 0 0 0] -> gives FE  coefficients without dt factor

function [bn] = ab3coefs(tn)


     tol=1e-7; %integration tolerance for quad

     %also should deal with degeneracies (FE)
     if tn(2) == tn(3) && tn(3)==tn(4)
         dt=tn(1)-tn(2);
         bn(1)=dt;
         bn(2)=0;
         bn(3)=0;
     elseif tn(2) ~=tn(3) && tn(3) == tn(4) %(AB2)
         %define lagrange interpolating polynomials
         l0=@(t) (t-tn(3))/(tn(2)-tn(3));
         l1=@(t) (t-tn(2))/(tn(3)-tn(2));

         %integrate inexactly with quad
         bn(1) = quad(l0,tn(2),tn(1),tol);
         bn(2) = quad(l1,tn(2),tn(1),tol);
         bn(3) = 0;
     else %general AB3 case:   
     %define lagrange interpolating polynomials
     l0=@(t) (t-tn(3))/(tn(2)-tn(3)).*(t-tn(4))/(tn(2)-tn(4));
     l1=@(t) (t-tn(2))/(tn(3)-tn(2)).*(t-tn(4))/(tn(3)-tn(4));
     l2=@(t) (t-tn(2))/(tn(4)-tn(2)).*(t-tn(3))/(tn(4)-tn(3));

     %integrate lagrange polynomials inexactly with quad
     bn(1) = quad(l0,tn(2),tn(1),tol);
     bn(2) = quad(l1,tn(2),tn(1),tol);
     bn(3) = quad(l2,tn(2),tn(1),tol);
     
end
