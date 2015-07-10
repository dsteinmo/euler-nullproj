%%
%
% Computes the Householder transformation of x into Hx.  That is, given a
% Householder vector v, this routine returns: 
%
%                       y = (I - 2*v*v')x. 
%
% y = apply_householder(x,v); 
function y = apply_householder(x,v)
    y = x - 2*v*(v'*x); 
end