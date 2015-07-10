function [rr,ss] = Nodes2Dquad(N)
    r = JacobiGL(0,0,N);
    [rr,ss] = meshgrid(r,r);
    rr=rr(:);
    ss=ss(:);
end