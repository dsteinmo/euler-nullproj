function pf2dquad(N,x,y,f)
    %note: the triangle interpolates to an 'order N' grid,
    %here I'm just using the order of the original function evaluation.
    
    %could generalize this by writing an InterpMatrix2Dquad function 
    %and using it here. For now, I'm lazy.
    global K Fmask
    
    %method='old';
    method='new';
    
    if strcmp(method,'old')
        

        %tic;
        hold on;
        for j=1:K
            xx = reshape(x(:,j),N+1,N+1);
            yy = reshape(y(:,j),N+1,N+1);
            ff = reshape(f(:,j),N+1,N+1);
            pcolor(xx,yy,ff);
        end
        hold off;
        shading flat;
        %toc;
    
    
    %new way: this totally wins for speed. hard-core. but,
    %only uses bilinear interpolation for within cells, so can look
    %"blocky".
    %suggest doing the interpolate to equispaced grid anyways.
    
    
    %figure;
    else
    %tic;
    %try constructing an 'fmask' that snakes its way thru all the nodes.

        Fmask2=Fmask;
        Fmask2(:,3)=flipud(Fmask2(:,3));
        Fmask2(:,4)=flipud(Fmask2(:,4));
        Fx= x(Fmask2(:),:);
        Fy= y(Fmask2(:),:);
        Ff= f(Fmask2(:),:);

        ph = patch(Fx,Fy,Ff);
        set(ph,'edgecolor','none');
    end
    axis tight;
    %colorbar;
    %toc;
    
end
