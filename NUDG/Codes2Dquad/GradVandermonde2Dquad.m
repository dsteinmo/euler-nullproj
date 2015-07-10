function [DVr,DVs] = GradVandermonde2Dquad(N,r,s)

% function [DVr,DVs] = GradVandermonde1D(N,r,s)
% Purpose : Initialize the gradient of the modal basis (i) at (r,s) at order N

DVr = zeros(length(r),(N+1)*(N+1));
DVs = zeros(length(s),(N+1)*(N+1));
% Initialize matrix
ind=1;
for i=0:N
    for j=0:N
       [DVr(:,ind)] = GradJacobiP(r(:),0,0,j).*JacobiP(s(:),0,0,i);
       [DVs(:,ind)] = JacobiP(r(:),0,0,j).*GradJacobiP(s(:),0,0,i);
       ind=ind+1;
    end
end
return
