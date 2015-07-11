function [uproj,vproj] = globalDivFreeProjection_curved(u,v,projMatrixStraight,rotDataStraight,projMatricesCurved)
    global curved straight
    [uproj(:,straight),vproj(:,straight)] = globalDivFreeProjection(u(:,straight),v(:,straight),projMatrixStraight,rotDataStraight);
    
    for n=1:length(curved)
        k1=curved(n);
        [uproj(:,k1),vproj(:,k1)] = divFreeProjection_curved(u(:,k1),v(:,k1),projMatricesCurved(:,:,n));
    end

end