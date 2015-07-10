function [mat,upsi,vpsi] = buildDivFreeMassMatrixquad(N,r,s,V,jac,rx,ry,sx,sy)
    %[psir,psis] = GradVandermonde2D(N-2,r,s);
    [psir,psis] = GradVandermonde2Dquad(N,r,s); %this seems like best bet.
    
    %[psir,psis] = GradVandermonde2D(N+1,r,s); %pretty shitty. -> dies in
                                               %real simulations
    %[psir,psis] = GradVandermonde2D(N+2,r,s); %definitely bad, yes - def!
    
    [Np, Npdiv] = size(psir);
    Npdiv = Npdiv-1;
    
    %Np = length(psir);
    %keyboard;
    
    if ~exist('jac','var') || ~exist('rx','var') ||  ~exist('ry','var') ...
            || ~exist('sx','var') || ~exist('sy','var')
        
        %if we're just on the reference triangle...
        jac=1;
        rx=ones(Np,1); ry=zeros(Np,1);
        sx=zeros(Np,1); sy=ones(Np,1);
    end
    
    %keyboard;   
    psix = repmat(rx,1,Npdiv+1).*psir + repmat(sx,1,Npdiv+1).*psis;
    psiy = repmat(ry,1,Npdiv+1).*psir + repmat(sy,1,Npdiv+1).*psis;
    
    upsi = -psiy;
    vpsi = psix;
    
    %upsi = -psis;
    %vpsi = psir; 
   
    upsi = upsi(:,2:end);
    vpsi = vpsi(:,2:end);

    %Npdiv = length(upsi)-1;

    mat = zeros(Npdiv,Npdiv);
    for ii=1:Npdiv
        for jj=1:Npdiv
            integrand = upsi(:,ii).*upsi(:,jj) + vpsi(:,ii).*vpsi(:,jj);
            %coefs = V\(integrand);
            %mat(ii,jj) = jac(1)*coefs(1)*(2/sqrt(2)); %integrate
            coefs = V\(jac.*integrand);
            mat(ii,jj) = coefs(1)*(2); %integrate %this line was only change for quads
        end
    end
end