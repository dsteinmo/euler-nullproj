%%
%
% Computes the Householder vector that transforms x --> ||x||e_j, where e_j
% is the standard Euclidean basis vector. 
%
%   v = get_householder(x,j);
%
% 27 Oct 2013 
% Sumedh Joshi
% Cornell University 

function v = get_householder(x,jj)
    e = 0*x; 
    e(jj) = 1.0; 
    v = x - norm(x)*e; 
    v = v/norm(v); 
end