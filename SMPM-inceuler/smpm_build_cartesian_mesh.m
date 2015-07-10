%%
%
% [x z] = smpm_build_cartesian_mesh(n,mx,mz,xlims,zlims,fpath,fname); 
%
%   Inputs
%       n     - number of GLL points per direction per subdomain. 
%       mx    - number of subdomains in the x-direction. 
%       mz    - number of subdomains in the z-direction.
%       xlims - [1x2]: bounds for the x-direction. 
%       zlims - [1x2]: bounds for the z-direction. 
%       fpath - the directory to save the file in. (optional). 
%       fname - the name to give the grid file.    (optional). 
%       
%
% Generate a cartesian (undeformed) mesh in xi-first indexing for reading
% by read_meshfile_data/header in the smpm code.  

%
% 7 Mar 2013
% Sumedh Joshi 

function [x z] = smpm_build_cartesian_mesh(n,mx,mz,xlims,zlims,fpath,fname)

    %
    % Get the GLL points. 
    xi  = lglnodes(n-1); 
    xi  = sort(xi); 
    eta = xi; 
    
    %
    % Build the 1D multidomain mesh in xi.
    Lz = (zlims(2) - zlims(1)) / mz; 
    xi = (xi + 1)*Lz/2 - (zlims(2)-zlims(1))/2;
    XI = repmat(xi,mz,1); 
    offset = repmat([0:mz-1]'*Lz,1,n)';
    offset = offset(:); 
    XI = XI + offset;
    
    %
    % Build the 1D multidomain mesh in eta. 
    Lx = (xlims(2) - xlims(1)) / mx; 
    eta = (eta + 1)*Lx/2 - (xlims(2)-xlims(1))/2; 
    ETA = repmat(eta,mx,1); 
    offset = repmat([0:mx-1]'*Lx,1,n)';
    offset = offset(:); 
    ETA = ETA + offset; 
    
    %
    % Do the meshgrid. 
    [x z] = meshgrid(ETA,XI); 
    x = x(:);
    z = z(:); 
    
    %
    % Add back in the offsets. 
    x = x + (xlims(2) - xlims(1))/2;
    z = z + (zlims(2) - zlims(1))/2;
    
    %
    % If asked, write to disk. 
    if nargin > 5
        smpm_write_meshfile(n,mx,mz,x,z,fpath,fname); 
    end
    
end