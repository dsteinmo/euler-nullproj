function BCType = CorrectBCTablequad(EToV,mapnodes,BCcode,K)

% function BCType = CorrectBCTable(EToV,BCType,mapnodes,BCcode);
% Purpose: Setup BCType for boundary conditions in 2D
% By Allan P. Engsig-Karup

%Globals2D;


%VNUM = [1 2;2 3;3 1]; % face orientations
VNUM = [1 2;2 3;3 4;4 1]; % face orientations

Nfaces=4;

BCType = zeros(K,Nfaces);

for k = 1:K    
    % Test for each edge
    for l = 1:Nfaces 
        m = EToV(k,VNUM(l,1)); n = EToV(k,VNUM(l,2));

        % if both points are on the boundary then it is a boundary point!
        ok=sum(ismember([m n],mapnodes));
        if ok==2 
            BCType(k,l)=BCcode; 
        end

    end
end
return