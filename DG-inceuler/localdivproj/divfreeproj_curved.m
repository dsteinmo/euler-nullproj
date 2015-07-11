%divFreeMassMat is (\vec{upsi}_i dot \vec{upsi}_j)

function [uproj,vproj] = divfreeproj_curved(u,v,upsi,vpsi,divFreeMassMat,cV,cW)
    

    sz = size(upsi);
    Npdiv = sz(2);
    
    rhs = zeros(Npdiv,1);
    for jj=1:Npdiv
        
        cu = cV*u;
        cv = cV*v;
        
        cupsi_j = cV*upsi(:,jj);
        cvpsi_j = cV*vpsi(:,jj);
        
        %integrand = u.*upsi(:,jj) + v.*vpsi(:,jj);
        
        cintegrand = cu.*cupsi_j + cv.*cvpsi_j;
        
        rhs(jj) = sum(cW.*cintegrand);
    end

    c = divFreeMassMat\rhs;
    uproj = 0*u;
    vproj = 0*v;
    for jj=1:Npdiv
        uproj = uproj+c(jj)*upsi(:,jj); 
        vproj = vproj+c(jj)*vpsi(:,jj);
    end
end