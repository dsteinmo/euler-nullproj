function pf2dquad_refel(N,x,y,f)
    %note: the triangle interpolates to an 'order N' grid,
    %here I'm just using the order of the original function evaluation.
    
    %could generalize this by writing an InterpMatrix2Dquad function 
    %and using it here. For now, I'm lazy.
    
    clf;
    
    xx = reshape(x(:),N+1,N+1);
    yy = reshape(y(:),N+1,N+1);
    ff = reshape(f(:),N+1,N+1);
    pcolor(xx,yy,ff); shading flat;

end