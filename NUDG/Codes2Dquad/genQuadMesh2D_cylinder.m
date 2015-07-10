%%generate simple uniform-element rectangules mesh on a rectangular domain
clear;
%%User-specified geometry details:
x1=0;
y1=.1;
Lx=5;  %domain length
Ly=.05;   %domain height (may be negative)

%specify number of elements in each direction
Nel_x = 40;
Nel_y = 10;
%%done user-specified info

%calculate number of elements
K=Nel_x*Nel_y;
%and opposite domain corner.
x2 = x1+Lx;
y2 = y1+Ly;

%fix this cross-referencing business
Nx=Nel_x;
Ny=Nel_y;

dx=Lx/Nx;
dy=Ly/Ny;

figure(1); clf;
Vnum=1;
Nv=(Nx+1)*(Ny+1);
VX=zeros(Nv,1);
VY=zeros(Nv,1);
%create list of vertices

boxleft = 1;
boxright= 1.5;
boxbot = 0.05;
boxtop = 0.10;
for j=0:Ny
    for i=0:Nx
        x=x1+(i)*dx;
        y=y1+(j)*dy;
        plot(x,y,'.');
        %text(x,y,num2str(Vnum));
        hold on;
        %pause;
        VX(Vnum)=x;
        VY(Vnum)=y;
        Vnum=Vnum+1;
    end
end

%create Element-To-Vertex connectivity table, we cheat by
%using the known structure of the vertices deployed above
Nfaces=4;
EToV = zeros(K,Nfaces); 
ind=1;
for i =1:Ny
    for j=1:Nx
        %do it with c.c.w orientation (note: this may break if you take x
        %or y negative!)
        EToV(ind,:) = (i-1)*(Nx+1)+ [j j+1 j+(Nx+2) j+(Nx+1)];
        plot(VX(EToV(ind,[1:Nfaces 1])),VY(EToV(ind,[1:Nfaces 1])),'-*r');
        hold on;
	ind=ind+1;
    end
end
%drawnow;
hold on;
%axis([0 5 0 0.15]);

%% Region 2 %%
Nv_old = Nv;
K_old = K;
Nx_old = Nx;
Ny_old = Ny;

x1=0;
y1=.05;
Lx=1;  %domain length
Ly=.05;   %domain height (may be negative)

%specify number of elements in each direction
Nel_x = 8;
Nel_y = 10;
%%done user-specified info

%calculate number of elements
K=Nel_x*Nel_y;
%and opposite domain corner.
x2 = x1+Lx;
y2 = y1+Ly;

%fix this cross-referencing business
Nx=Nel_x;
Ny=Nel_y;

dx=Lx/Nx;
dy=Ly/Ny;

figure(1);
Nv=(Nx+1)*(Ny+1);

for j=0:Ny-1  %stop at Ny-1 to skip the last row, it's already in the mesh.
    for i=0:Nx
        x=x1+(i)*dx;
        y=y1+(j)*dy;
        plot(x,y,'.');
        %text(x,y,num2str(Vnum));
        hold on;
        %pause;
        VX(Vnum)=x;
        VY(Vnum)=y;
        Vnum=Vnum+1;
    end
end

%create Element-To-Vertex connectivity table, we cheat by
%using the known structure of the vertices deployed above
Nfaces=4;
EToV = [EToV; zeros(K,Nfaces)];
%VX = [VX; zeros(Nv,1)];
%VY = [VY; zeros(Nv,1)];


%ind=K_old+1;
for i =1:Ny-1 %skip top row.
    for j=1:Nx
        %do it with c.c.w orientation (note: this may break if you take x
        %or y negative!)
        EToV(ind,:) = Nv_old+(i-1)*(Nx+1)+ [j j+1 j+(Nx+2) j+(Nx+1)];
        plot(VX(EToV(ind,[1:Nfaces 1])),VY(EToV(ind,[1:Nfaces 1])),'-*r');
        %hold on;
	    ind=ind+1;
    end
end

%Stitch to Region 1 for the top row
for j=1:Nx
    EToV(ind,:) = [Nv_old + (Ny-1)*(Nx+1)+j Nv_old + (Ny-1)*(Nx+1)+j+1 1+j j];
    plot(VX(EToV(ind,[1:Nfaces 1])),VY(EToV(ind,[1:Nfaces 1])),'-*r');
    hold on;
    ind=ind+1;
end

%drawnow;
%axis([0 5 0 0.15]);
%axis([0 2 0.05 0.15])

%% Region 3 %%
Nv_old = Nv + Nv_old;
K_old = K + K_old;
Nx_old = Nx+Nx_old;
Ny_old = Ny+Ny_old;

x1=2;
y1=.05;
Lx=3;  %domain length
Ly=.05;   %domain height (may be negative)

%specify number of elements in each direction
Nel_x = 8*3;
Nel_y = 10;
%%done user-specified info

%calculate number of elements
K=Nel_x*Nel_y;
%and opposite domain corner.
x2 = x1+Lx;
y2 = y1+Ly;

%fix this cross-referencing business
Nx=Nel_x;
Ny=Nel_y;

dx=Lx/Nx;
dy=Ly/Ny;

figure(1);
Nv=(Nx+1)*(Ny+1);

for j=0:Ny-1  %stop at Ny-1 to skip the last row, it's already in the mesh.
    for i=0:Nx
        x=x1+(i)*dx;
        y=y1+(j)*dy;
        plot(x,y,'.');
        %text(x,y,num2str(Vnum));
        hold on;
        %pause;
        VX(Vnum)=x;
        VY(Vnum)=y;
        Vnum=Vnum+1;
    end
end

%create Element-To-Vertex connectivity table, we cheat by
%using the known structure of the vertices deployed above
Nfaces=4;

for i =1:Ny-1 %skip top row.
    for j=1:Nx
        %do it with c.c.w orientation (note: this may break if you take x
        %or y negative!)
        EToV(ind,:) = (Nv_old)+(i-1)*(Nx+1)-9+ [j j+1 j+(Nx+2) j+(Nx+1)];
        plot(VX(EToV(ind,[1:Nfaces 1])),VY(EToV(ind,[1:Nfaces 1])),'-*r');
        hold on;
	    ind=ind+1;
    end
end

%Stitch to Region 1 for the top row
for j=1:Nx
    EToV(ind,:) = [Nv_old + (Ny-1)*(Nx+1)+j-9 Nv_old + (Ny-1)*(Nx+1)+j+1-9 1+j+16 j+16];
    plot(VX(EToV(ind,[1:Nfaces 1])),VY(EToV(ind,[1:Nfaces 1])),'-*r');
    hold on;
    ind=ind+1;
end

%drawnow;
%axis([0 5 0.05 0.15])

%%Region 4%%
Nv_old = Nv + Nv_old;
K_old = K + K_old;
Nx_old = Nx+Nx_old;
Ny_old = Ny+Ny_old;

x1=0;
y1=0;
Lx=5;  %domain length
Ly=.05;   %domain height (may be negative)

%specify number of elements in each direction
Nel_x = 40;
Nel_y = 10;
%%done user-specified info

%calculate number of elements
K=Nel_x*Nel_y;
%and opposite domain corner.
x2 = x1+Lx;
y2 = y1+Ly;

%fix this cross-referencing business
Nx=Nel_x;
Ny=Nel_y;

dx=Lx/Nx;
dy=Ly/Ny;

figure(1);
Nv=(Nx+1)*(Ny+1);

for j=0:Ny-1  %stop at Ny-1 to skip the last row, it's already in the mesh.
    for i=0:Nx
        x=x1+(i)*dx;
        y=y1+(j)*dy;
        plot(x,y,'.');
        %text(x,y,num2str(Vnum));
        hold on;
        %pause;
        VX(Vnum)=x;
        VY(Vnum)=y;
        Vnum=Vnum+1;
    end
end

for i =1:Ny-1 %skip top row.
    for j=1:Nx
        %do it with c.c.w orientation (note: this may break if you take x
        %or y negative!)
        EToV(ind,:) = (Nv_old)+(i-1)*(Nx+1)-34 + [j j+1 j+(Nx+2) j+(Nx+1)];
        plot(VX(EToV(ind,[1:Nfaces 1])),VY(EToV(ind,[1:Nfaces 1])),'-*r');
        hold on;
	    ind=ind+1;
    end
end



%Stitch to Region 1 for the top row
for j=1:8
    EToV(ind,:) = [Nv_old + (Ny-1)*(Nx+1)+j-34 Nv_old + (Ny-1)*(Nx+1)+j-33   1+j+451 j+451];
    plot(VX(EToV(ind,[1:Nfaces 1])),VY(EToV(ind,[1:Nfaces 1])),'-*r');
    hold on;
    ind=ind+1;
end

%Stitch to Region 1 for the top row
for j=17:40
    EToV(ind,:) = [Nv_old + (Ny-1)*(Nx+1)+j-34 Nv_old + (Ny-1)*(Nx+1)+j-33   525+j+1 j+525];
    plot(VX(EToV(ind,[1:Nfaces 1])),VY(EToV(ind,[1:Nfaces 1])),'-*r');
    hold on;
    ind=ind+1;
end

axis([0 5 0 0.15])
drawnow;

bottomBdry = (find(abs(VY-0)<1e-4))';

bottomEdge = [(1:length(bottomBdry)-1)' (2:length(bottomBdry))']; 

rightBdry = [832 873 914 955 996 1037 1078 1119 1160 1201 566 591 616 641 666 691 716 741 766 791 41   82   123  164  205  246 287 328 369 451];

rightEdge = length(bottomEdge)+[(1:length(rightBdry)-1)' (2:length(rightBdry))'];

topBdry = 451:-1:411;

topEdge = length(bottomEdge)+length(rightEdge)+[(1:length(topBdry)-1)' (2:length(topBdry))'];

leftBdry = [411 370 329 288 247 206 165 124 83 42 1 533 524 515 506 497 488 479 470 461 452 1161 1120 1079 1038 997 956 915 874 833 792];

leftEdge = length(bottomEdge)+length(rightEdge)+length(topEdge)+[(1:length(leftBdry)-1)' (2:length(leftBdry))'];

cylBdry = [9 10 11 12 13 14 15 16 17 767 742 717 692 667 642 617 592 567 542 1177 1176 1176 1174 1173 1172 1171 1170 1169 460 469 478 487 496 505 514 523 532 541 9];

cylEdge = length(bottomEdge)+length(rightEdge)+length(topEdge)+length(leftEdge)+[(1:length(cylBdry)-1)' (2:length(cylBdry))'];

bdryNodes = [bottomBdry rightBdry topBdry leftBdry cylBdry];

%bdryEdge = [bottomEdge; rightEdge; topEdge; leftEdge; cylEdge];

bdryEdge = [(1:length(bdryNodes)-1)' (2:length(bdryNodes))'];

%bdryNodes = [(1:Nx) ((Nx+1):(Nx+1):Nv) (Nv-1):-1:(Nv-Nx) (Nv-Nx-(Nx+1)):-(Nx+1):1];
%bdryNodes = [(1:Nx) ((Nx+1):(Nx+1):Nv) (Nv-Nx-(Nx+1)):-(Nx+1):1 (Nv-Nx):Nv+(Nx/5)];

%edge connectivity into these nodes is most trivial:
%bdryEdge = [(1:length(bdryNodes)-1)' (2:length(bdryNodes))'];

hold on;

%plot(VX(bdryNodes),VY(bdryNodes),'-k','linewidth',2);

%plot(VX(bdryNodes(bdryEdge(:))),VY(bdryNodes(bdryEdge(:))),'-k','linewidth',2); 
