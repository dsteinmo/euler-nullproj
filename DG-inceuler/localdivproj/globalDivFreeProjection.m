function [uproj,vproj] = globalDivFreeProjection(u,v,projectionMatrix,rotData)
    %rotate (u,v) to (ur,us) coordinates
    
    %old code:
    %ur = rotData.Ainv11.*(u-rotData.x0) + rotData.Ainv12.*(v-rotData.y0); 
    %us = rotData.Ainv21.*(u-rotData.x0) + rotData.Ainv22.*(v-rotData.y0); 
    %project on standard element (i.e., call local projection)
    %[urproj,usproj] = divFreeProjection(ur,us,projectionMatrix);
    %rotate back
    %uproj = rotData.A11.*urproj + rotData.A12.*usproj + rotData.x0;
    %vproj = rotData.A21.*urproj + rotData.A22.*usproj + rotData.y0;
    
    %new code: suggest doing this without the shift:
    %rotate (u,v) to (ur,us) coordinates
    ur = rotData.Ainv11.*(u) + rotData.Ainv12.*(v); 
    us = rotData.Ainv21.*(u) + rotData.Ainv22.*(v); 
    %project on standard element (i.e., call local projection)
    [urproj,usproj] = divFreeProjection(ur,us,projectionMatrix);
    %rotate back
    uproj = rotData.A11.*urproj + rotData.A12.*usproj;
    vproj = rotData.A21.*urproj + rotData.A22.*usproj;
end