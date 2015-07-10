%purpose: integrate a field 'f' over an element or over the domain
%uses a cubature structure 'cub' that can be built by calling
%'cub = CubatureVolumeMesh2D(CubatureOrder);'

%all three if-else cases have been validated for 
%integrand =1 case.
function result = dgint2d_curved(f,cub,k1)
    %interpolate function to cubature nodes
    cf = cub.V*f;
    
    %single-element case
    if exist('k1','var')
        if strcmp(k1,'ref') %reference element case.
            result = sum((cub.w).*cf);
        else
            result = sum(cub.W(:,k1).*cf); %note cub.W contains jacobian
        end
    else
        %whole-domain case        
        result = sum(sum((cub.W).*cf)); %note cub.W contains jacobian
    end
end