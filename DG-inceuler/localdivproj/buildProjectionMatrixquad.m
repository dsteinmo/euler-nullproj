function projectionMatrix = buildProjectionMatrixquad(N,r,s,V)

    %To retain stability, need to take N_psi = N here.
    N_psi = N;
   
    Np = length(r);
    [psir,psis] = GradVandermonde2Dquad(N_psi,r,s);
    upsi = -psis(:,2:end);
    vpsi = psir(:,2:end);
    %basis is genuinely divergence-free. (I checked)
    
    
    divFreeMassMat = buildDivFreeMassMatrixquad(N_psi,r,s,V);

    projectionMatrix = zeros(2*Np,2*Np);
    g = zeros(2*Np,1);
    for ii=1:2*Np;
        g(ii)=1;
        u = g(1:end/2);
        v = g(end/2+1:end);
        [uproj,vproj] = divfreeprojquad(u,v,upsi,vpsi,divFreeMassMat,V);

        col = [uproj;vproj];
        projectionMatrix(:,ii) = col;

        g(ii)=0;
    end
end