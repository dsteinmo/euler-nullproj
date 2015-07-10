

function [mat,upsi,vpsi] = buildDivFreeMassMatrix_curved(Dx,Dy,V,cV,cW)
       %V = psi
       upsi = -Dy*V;
       vpsi = Dx*V;
       %can use my 'monomials' function if this turns out to be bad.
       %[upsi,vpsi] = buildDivFreeMonomials(...)
       Np = length(V);
       
       upsi = upsi(:,2:end);
       vpsi = vpsi(:,2:end);
       
       Npdiv = Np-1;
       
       mat = zeros(Npdiv,Npdiv);
       for ii=1:Npdiv
           for jj=1:Npdiv
               %interpolate to cubature nodes
               cupsi_i = cV*upsi(:,ii);
               cupsi_j = cV*upsi(:,jj);
               
               cvpsi_i = cV*vpsi(:,ii);
               cvpsi_j = cV*vpsi(:,jj);
               
               cintegrand = cupsi_i.*cupsi_j + cvpsi_i.*cvpsi_j;
               
               %original thought:
               %integrand = upsi(:,ii).*upsi(:,jj) + vpsi(:,ii).*vpsi(:,jj);
               %cintegrand = cV*integrand;
               
               
               %integrate
               mat(ii,jj) = sum(cW.*cintegrand);
           end
       end
       
end