%%generate simple uniform-element rectangules mesh on a rectangular domain
clear;
%%User-specified geometry details:
x1=0;
y1=0;
Lx=5;  %domain length
Ly=.15;   %domain height (may be negative)

%specify number of elements in each direction
Nel_x = 10;
Nel_y = 200;
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
for j=0:Ny
    for i=0:Nx
        x=x1+(i)*dx;
        y=y1+(j)*dy;
%         plot(x,y,'.');
%         text(x,y,num2str(Vnum));
%         hold on;
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
for i =1:Ny
    for j=1:Nx
        %do it with c.c.w orientation (note: this may break if you take x
        %or y negative!)
        EToV(i,:) = (i-1)*(Nx+1)+ [j j+1 j+(Nx+2) j+(Nx+1)];
        plot(VX(EToV(i,[1:Nfaces 1])),VY(EToV(i,[1:Nfaces 1])),'-*r');
        hold on;
    end
end
drawnow;
hold off;

bdryNodes = [(1:Nx) ((Nx+1):(Nx+1):Nv) (Nv-1):-1:(Nv-Nx) (Nv-Nx-(Nx+1)):-(Nx+1):1];
%edge connectivity into these nodes is most trivial:
bdryEdge = [(1:length(bdryNodes)-1)' (2:length(bdryNodes))'];

%figure(2)
%plot(VX(bdryNodes),VY(bdryNodes),'.-g')
% x=linspace(x1,x2,Nel_x+1);
% y=linspace(y1,y2,Nel_y+1);
% 
% [xx,yy]=meshgrid(x,y);
% 
% %make global list of vertices
% VX= xx(:);
% VY= yy(:);
% 
% 
% plot(VX,VY,'.');
% hold on;
% for j=1:length(VX)
%     plot(VX(j),VY(j),'*r');
%     pause;
% end