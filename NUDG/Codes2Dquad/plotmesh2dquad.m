function plotmesh2dquad()
    global VX VY EToV
    
    hold on;
    xmin = min(VX);
    xmax = max(VX);
    ymin = min(VY);
    ymax = max(VY);
    Lx = xmax-xmin;
    Ly = ymax-ymin;
    axis([xmin-0.1*Lx xmax+0.1*Lx ymin-0.1*Ly ymax+0.1*Ly]);
    for j=1:length(EToV)
        plot(VX([EToV(j,:) EToV(j,1)]),VY([EToV(j,:) EToV(j,1)]),'-k');
    end
    hold off; 
end