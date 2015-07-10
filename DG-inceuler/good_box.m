clear;

node = [0 0; 
        1 0;
        1 1;
        0 1];
        

edge = [1 2;
        2 3;
        3 4;
        4 1];

hdata.hmax = .05;
    
[Vert,EToV] = mesh2d(node,edge,hdata);

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
    
    
    %refinement stuff.
    refineflag =1;
    while refineflag == 1
        disp('Zoom to region that needs refinement. Press any key when done.');
        pause;
        
        disp('select an element for refinement (click vertex opposite longest edge)...');
        [myx,myy] = ginput(1);
        mytol=.05; %10m tolerance
        dists = sqrt((VX-myx).^2 + (VY-myy).^2);
        myV = find(dists<mytol);
        %myV
        
        if length(myV) ~= 1
            disp('No vertex selected, assuming you''re done refining...');
            refineflag=0;
        else
            disp(['You clicked on vertex number: ' num2str(myV)]);
            myelement = mod(find(EToV==myV),K);
            if myelement == 0
                myelement = K;
            end
            if length(myelement) ~= 1
                disp('Error: You picked a vertex with multiple elements. Click corner vertices only. (for now)');
            else
                disp(['Element number: ' num2str(myelement)]);    
                localEToV = EToV(myelement,:); %1x3 array of the element we want to refine
                myVind = find(localEToV==myV); %index of myV in localEToV
                otherinds = find(localEToV ~= myV); %other two vertex indices
                other(1) = localEToV(otherinds(1)); %get global vertex numbers
                other(2) = localEToV(otherinds(2)); %of other two

                %calculate midpoint of other two vertices
                midpointx = 0.5*(VX(other(1)) + VX(other(2)));
                midpointy = 0.5*(VY(other(1)) + VY(other(2)));

                %add midpoint as new vertex to global vertex table
                VX=[VX;midpointx]; VY=[VY;midpointy];
                %midpoint's is at the end of the table, so its index is
                %length(table)
                midind = length(VX);

                %draw new line segment on mesh to indicate that refinement is
                %being done
                hold on;
                plot([midpointx VX(myV)],[midpointy VY(myV)],'-b'); drawnow;
                hold off;

                %find the element our 'split-face' is shared with and refine it
                %appropriately
                myinds = mod(find(EToV == other(1)),K);
                myinds(myinds==0) = K; %need this hack b/c of matlab.
                %if myinds == 0 
                %    myinds =K;
                %end
                
                myinds = myinds(myinds ~= myelement); %remove the element we started with from the contenders
                numpossible = length(myinds);
                possibleelements = EToV(myinds,:);
                theind = mod(find(possibleelements == other(2)),numpossible);
                if theind == 0 %need this hack to get around matlab no zero index problems
                    theind = numpossible;
                end
                otherelement = myinds(theind);
                othervertices = EToV (otherelement,:);       %global vertex numbers of the 'other element'
                %othervertices = possibleelements(theind,:); %global vertex numbers of the 'other element'

                tmpinds = find(othervertices ~= other(2));
                newmyVind = othervertices(tmpinds) ~= other(1);
                newmyVind = tmpinds(newmyVind);

                newmyV = othervertices(newmyVind);

                %draw new line segment on mesh to indicate that refinement is
                %being done
                hold on;
                plot([midpointx VX(newmyV)],[midpointy VY(newmyV)],'-b'); drawnow;
                hold off;

                %modify old element indices to include midpoint
                EToV(myelement,:) = [midind other(1) myV];
                %add a new element to the connectivity table
                newel = [midind myV other(2)];
                EToV = [EToV; newel];

                %modify old element indices to include midpoint
                EToV(otherelement,:) = [newmyV other(1) midind];
                %add a new element to the connectivity table
                newel = [newmyV midind other(2)];
                EToV = [EToV; newel];


                %increment number of vertices and increment twice the of elements
                Nv=Nv+1
                K=K+2

                %highlight our refined elemenets (repeat first vertex to close
                %off triangles)
                hold on;
                plot(VX([EToV(myelement,:) EToV(myelement,1)]),VY([EToV(myelement,:) EToV(myelement,1)]),'-r','linewidth',2);
                plot(VX([EToV(otherelement,:) EToV(otherelement,1)]),VY([EToV(otherelement,:) EToV(otherelement,1)]),'-g','linewidth',2);
                plot(VX([EToV(end-1,:) EToV(end-1,1)]),VY([EToV(end-1,:) EToV(end-1,1)]),'-m','linewidth',2);
                plot(VX([EToV(end,:) EToV(end,1)]),VY([EToV(end,:) EToV(end,1)]),'-c','linewidth',2);
                drawnow;
                hold off;

                disp('done refining that element and its affected neighbour.');
              

            end

        end 
    end %end refinement stuff while loop
    
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
    if exist('N','var')
       dof = nchoosek(N+2,2)*K;
       disp(['DG simulation with order N=' num2str(N) ' basis polynomials will have ' num2str(dof) ' degrees of freedom.']);
    end
    
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
    BCType = CorrectBCTable(EToV,BCType,nodesOuter,Wall,K);

    %Need to do this to make vertex arrays consistent with main scripts.
    VX = VX';
    VY = VY';
