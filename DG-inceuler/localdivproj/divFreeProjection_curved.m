%this is to be called for single elements only,
%each with their own uniquely-defined projection matrix.
function [uproj,vproj] = divFreeProjection_curved(u,v,projectionMatrix)
    invec = [u;v];
    outvec = projectionMatrix*invec;
    uproj = outvec(1:end/2,:); vproj = outvec(end/2+1:end,:);
end