function projectionMatrix = buildProjectionMatrix_curved(Dx,Dy,V,cV,cW)

    Np = length(V);
        
    [divFreeMassMat,upsi,vpsi] = buildDivFreeMassMatrix_curved(Dx,Dy,V,cV,cW);

    projectionMatrix = zeros(2*Np,2*Np);
    
    %build local operator the cheap way by passing in unit vectors.
    g = zeros(2*Np,1);
    for ii=1:2*Np;
        g(ii)=1;
        u = g(1:end/2);
        v = g(end/2+1:end);
        [uproj,vproj] = divfreeproj_curved(u,v,upsi,vpsi,divFreeMassMat,cV,cW);

        col = [uproj;vproj];
        projectionMatrix(:,ii) = col;

        g(ii)=0;
    end
end