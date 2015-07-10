%build coordinate rotation information on all elements.
%for going from (x,y) to (r,s) coordinates

%returns as a matlab struct
function r = buildCoordRotationData
    global K EToV VX VY Np

    r.A11 = zeros(Np,K); r.A12 = zeros(Np,K);
    r.A21 = zeros(Np,K); r.A22 = zeros(Np,K);
    r.Ainv11 = zeros(Np,K); r.Ainv12 = zeros(Np,K);
    r.Ainv21 = zeros(Np,K); r.Ainv22 = zeros(Np,K);
    r.x0 = zeros(Np,K); r.y0 = zeros(Np,K);
    %ten extra fields to store, altogether.
    for k1=1:K
        va = EToV(k1,1)'; vb = EToV(k1,2)'; vc = EToV(k1,3)';
        xa = VX(va); xb = VX(vb); xc = VX(vc);
        ya = VY(va); yb = VY(vb); yc = VY(vc);

        %x = 0.5*(-(r+s)*xa+(1+r)*xb+(1+s)*xc);
        %y = 0.5*(-(r+s)*ya+(1+r)*yb+(1+s)*yc);
        %or, 
        %(x,y) = A *(r,s) + 0.5*(xb+xc,yb+yc)
        %where:
        Amat = 0.5*[(xb-xa) (xc-xa);
                    (yb-ya) (yc-ya)];
                
        Ainv = inv(Amat);
        
        r.A11(:,k1) = repmat(Amat(1,1),Np,1); r.A12(:,k1) = repmat(Amat(1,2),Np,1);
        r.A21(:,k1) = repmat(Amat(2,1),Np,1); r.A22(:,k1) = repmat(Amat(2,2),Np,1);
        r.Ainv11(:,k1) = repmat(Ainv(1,1),Np,1); r.Ainv12(:,k1) = repmat(Ainv(1,2),Np,1);
        r.Ainv21(:,k1) = repmat(Ainv(2,1),Np,1); r.Ainv22(:,k1) = repmat(Ainv(2,2),Np,1);

        %may not need these two given recent changes.
        r.x0(:,k1) = repmat(0.5*(xb+xc),Np,1);
        r.y0(:,k1) = repmat(0.5*(yb+yc),Np,1);
        
    end
end