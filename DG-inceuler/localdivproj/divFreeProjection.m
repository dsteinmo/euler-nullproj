%this is meant for use on the reference triangle only
%so it is called by globalDivFreeProjection on (u,v)
%once rotated to the (r,s) coordinate frame.
function [uproj,vproj] = divFreeProjection(u,v,projectionMatrix)
    invec = [u;v];
    outvec = projectionMatrix*invec;
    uproj = outvec(1:end/2,:); vproj = outvec(end/2+1:end,:);
end