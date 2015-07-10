function contour2dquad(N,x,y,f,conts)
    %note: the triangle interpolates to an 'order N' grid,
    %here I'm just using the order of the original function evaluation.
    
    %could generalize this by writing an InterpMatrix2Dquad function 
    %and using it here. For now, I'm lazy.
    global K Fmask
    
    hold on;
    for j=1:K
        xx = reshape(x(:,j),N+1,N+1);
        yy = reshape(y(:,j),N+1,N+1);
        ff = reshape(f(:,j),N+1,N+1);
        contour(xx,yy,ff,conts,'-k');
    end
    hold off;
    axis tight;     
end
