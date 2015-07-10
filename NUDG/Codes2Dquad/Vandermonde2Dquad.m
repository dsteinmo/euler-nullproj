function [V2D] = Vandermonde2Dquad(N,r,s)

% function [V2D] = Vandermonde2Dquad(N,r)
% Purpose : Initialize the 2D Vandermonde Matrix, V_{ij} = phi_i(r_k) phi_j(s_k);

V2D = zeros(length(r),(N+1)^2);
ind=1;
for i=1:N+1
    for j=1:N+1
        V2D(:,ind) = JacobiP(s(:), 0, 0, i-1).*JacobiP(r(:), 0, 0, j-1);
        ind=ind+1;
    end
end
%Q: do these need to be scaled somehow?
return
