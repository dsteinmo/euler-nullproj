%analytically differentiations x^i*y^j basis functions
%for given (x,y) coords such that (u,v) is a non-divergent basis
%N is polynomial order used in DG simulations 

%meant for a single element.
%this assumes triangle, so Np = n+2 choose 2 and inner loop is only up to
%i

function [u,v] = buildDivFreeMonomials(N,x,y)
    %Np = nchoosek(N+2,2);
    Np = length(x);
    n=1;
    u = zeros(Np,Np);
    v = zeros(Np,Np);
    for i=0:N
        for j=0:i
            u(:,n) = -(j*x.^(i-j).*(y.^(j-1)));
            v(:,n) = (i-j)*x.^(i-j-1).*(y.^j);
            n=n+1;
        end
    end
        
end