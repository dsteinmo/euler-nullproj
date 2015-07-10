clear;

%(x-5/a)^2 + (y-.5/b)^2 = 1;
    

% edge = [1 2;
%         2 3;
%         3 4;
%         4 1];

%hdata.hmax = .15;
hdata.hmax = .12;
Nt=100;
dt=2*pi/Nt;
t=(0:Nt-1)*dt;
a=5; %major semi-axis
b=.5;
x0=5;
y0=-.5;
x=x0+a*cos(t);
y=y0+b*sin(t);
node=[x' y'];
edge = [(1:length(node))' [2:length(node) 1]'];
[Vert,EToV] = mesh2d(node,edge,hdata);
%[Vert,EToV] = mesh2d(node,[],hdata);


    axis on;
    grid on; 
    xlabel('x (m)'); ylabel('y (m)');
    axis tight;
    dimp = size(Vert); dimt = size(EToV);
    Nv = dimp(1); K = dimt(1);
    
    
    disp(['Mesh generated. ' num2str(Nv) ' vertices on ' num2str(K) ' elements.']);
    
    
    %stuff below for DG
    VX = Vert(:,1); VY = Vert(:,2);
    
    % Reorder elements to ensure counter clockwise orientation
    ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
    bx = VX(EToV(:,2)); by = VY(EToV(:,2));
    cx = VX(EToV(:,3)); cy = VY(EToV(:,3));

    D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
    i = find(D<0);
    EToV(i,:) = EToV(i,[1 3 2]);
    %done reordering

    % Build connectivity matrix
    %[EToE, EToF] = tiConnect2D(EToV);
    
    
    %axis([1.1975e4-500 1.1975e4+500 1.1853e4-500 1.1853e4+500]);
    
    %axis([7800 8700 4300 5300]);
 
    
    % Reorder elements to ensure counter clockwise orientation
    ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
    bx = VX(EToV(:,2)); by = VY(EToV(:,2));
    cx = VX(EToV(:,3)); cy = VY(EToV(:,3));

    D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
    i = find(D<0);
    EToV(i,:) = EToV(i,[1 3 2]);
    %done reordering

    % Build connectivity matrix
    [EToE, EToF] = tiConnect2D(EToV);
    
    %calculate # of DoF's if we were asked to
%     if exist('N','var')
%        dof = nchoosek(N+2,2)*K;
%        disp(['DG simulation with order N=' num2str(N) ' basis polynomials will have ' num2str(dof) ' degrees of freedom.']);
%     end
    
    %Now set up boundary tables.
    %find all boundary nodes (nodesOuter). Ain't nothin to it, but to do it.
    tol = 1e-8;
    edgenum = findedge(Vert,node,edge,tol); %mesh2d routine
    kk=1;
    for jj = 1:length(edgenum)
        %if edgenum of point jj is nonzero, then it lies on the boundary,
        %so put it in list of boundary pts.
        if edgenum(jj) ~= 0  
            nodesOuter(kk) = jj;
            kk=kk+1;
        end        
    end
    
    %allocate BCType table.
    BCType = 0*EToE;
    
    %Insert the correct BC codes for boundaries
    Wall=3;
    BCType = CorrectBCTable_derek(EToV,VX,VY,node,edge,BCType,Wall);


    %Need to do this to make vertex arrays consistent with main scripts.
    VX = VX';
    VY = VY';
