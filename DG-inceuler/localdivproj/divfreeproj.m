%divFreeMassMat is (\vec{upsi}_i dot \vec{upsi}_j)

function [uproj,vproj] = divfreeproj(u,v,upsi,vpsi,divFreeMassMat,V,jac)
    if ~exist('jac','var')
        jac=1;
    end

    sz = size(upsi);
    Npdiv = sz(2);
    
    rhs = zeros(Npdiv,1);
    for jj=1:Npdiv
        integrand = u.*upsi(:,jj) + v.*vpsi(:,jj);
        %coefs = V\integrand;
        %rhs(jj) = jac(1)*coefs(1)*(2/sqrt(2));
        
        coefs = V\(jac.*integrand);
        rhs(jj) = coefs(1)*(2/sqrt(2));
    end

    c = divFreeMassMat\rhs;
    uproj = 0*u;
    vproj = 0*v;
    for jj=1:Npdiv
        uproj = uproj+c(jj)*upsi(:,jj); 
        vproj = vproj+c(jj)*vpsi(:,jj);
    end
end