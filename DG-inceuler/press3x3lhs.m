%for implicit time-stepping
%impl_param = 1   => backward euler
%impl_param =.5   => trapezoid rule (if RHS is appropriately modified)
%impl_param =2/3  => BDF2 (if RHS is appropriately modified)
%impl_param =6/11 => BDF3 (if RHS is appropriately modified)

function result = press3x3lhs(invec,dt,c0,impl_param)

    if ~exist('impl_param','var')
        impl_param=1; %for BE step.
    end
    
    Globals2D;
    u=invec(1:Np*K); u=reshape(u,Np,K);
    v=invec(Np*K+1:2*Np*K); v=reshape(v,Np,K);
    p=invec(2*Np*K+1:3*Np*K); p=reshape(p,Np,K);
    
    [px,py]=Grad2D(p);
    

    du = zeros(3*Nfp,K);
    dv = zeros(3*Nfp,K);
    dp = zeros(3*Nfp,K);
    
    uM=u(vmapM); uP=u(vmapP);
    vM=v(vmapM); vP=v(vmapP);
    pM=p(vmapM); pP=p(vmapP);
    
    %ghost the no normal flow condition at Walls
    uP(mapW) = uM(mapW) - 2*nx(mapW).*(nx(mapW).*uM(mapW) + ny(mapW).*vM(mapW));
    vP(mapW) = vM(mapW) - 2*ny(mapW).*(nx(mapW).*uM(mapW) + ny(mapW).*vM(mapW));
    

    lambda=c0; %assumes constant, which is should be
    c02=c0^2;
    
    %gtau=8e4;
    %gtau=1e-2;
    %gtau=8e-4;
    %gtau=100;
    %gtau=100;
    %gtau=10;
    
    %gtau=0;
    
    %try a fully SIPDG inspired method
    % Compute minimum height of elements either side of each edge
    hmin = min(2*J(vmapP)./sJ(mapP), 2*J(vmapM)./sJ(mapM));
    tau = reshape(Np./hmin, Nfp*Nfaces, K); 
    
    du(:)=uM-uP;
    dv(:)=vM-vP;
    dp(:)=pM-pP;
    
    %lhs1=u+impl_param*dt*(px-LIFT*(Fscale.*(nx.*dp-lambda.*du)/2)); %old
    %lhs2=v+impl_param*dt*(py-LIFT*(Fscale.*(ny.*dp-lambda.*dv)/2)); %old
    %lhs1=u+impl_param*dt*(px-LIFT*(Fscale.*(nx.*dp-gtau.*du)/2)); 
    %lhs2=v+impl_param*dt*(py-LIFT*(Fscale.*(ny.*dp-gtau.*dv)/2)); 
    %artificial compressibilitiy
    %keyboard;
    %lhs3=p+impl_param*dt*(c02*Div2D(u,v)-LIFT*(Fscale.*(c02*nx.*du+c02*ny.*dv -lambda.*dp)/2));
    %lhs3= (c02*Div2D(u,v)-LIFT*(Fscale.*(c02*nx.*du+c02*ny.*dv -lambda.*dp)/2)); %old
    %lhs3= (Div2D(u,v)-LIFT*(Fscale.*(nx.*du+ny.*dv -(1/gtau).*dp)/2)); %new
    
    lhs1=u+impl_param*dt*(px-LIFT*(Fscale.*(nx.*dp)/2)); 
    lhs2=v+impl_param*dt*(py-LIFT*(Fscale.*(ny.*dp)/2));
    lhs3= (Div2D(u,v)-LIFT*(Fscale.*(nx.*du+ny.*dv +tau.*dp)/2)); %new
    
    result=[lhs1(:);lhs2(:);lhs3(:)];
     
end
