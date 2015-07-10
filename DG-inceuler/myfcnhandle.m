
function result = myfcnhandle(invec,dt,c0)
    Globals2D;
    u=invec(1:Np*K); u=reshape(u,Np,K);
    v=invec(Np*K+1:2*Np*K); v=reshape(v,Np,K);
    p=invec(2*Np*K+1:3*Np*K); p=reshape(p,Np,K);
    
    [px,py]=Grad2D(p);
    

    du = zeros(3*Nfp,K);
    dv = zeros(3*Nfp,K);
    dp=zeros(3*Nfp,K);
    
    uM=u(vmapM); uP=u(vmapP);
    vM=v(vmapM); vP=v(vmapP);
    pM=p(vmapM); pP=p(vmapP);
    
    %ghost the no normal flow condition at Walls
    uP(mapW) = uM(mapW) - 2*nx(mapW).*(nx(mapW).*uM(mapW) + ny(mapW).*vM(mapW));
    vP(mapW) = vM(mapW) - 2*ny(mapW).*(nx(mapW).*uM(mapW) + ny(mapW).*vM(mapW));
    

    lambda=c0; %assumes constant, which is should be
    c02=c0^2;
    
    du(:)=uM-uP;
    dv(:)=vM-vP;
    dp(:)=pM-pP;
    
    lhs1=u+dt*(px-LIFT*(Fscale.*(nx.*dp-lambda.*du)/2));
    lhs2=v+dt*(py-LIFT*(Fscale.*(ny.*dp-lambda.*dv)/2));
    %artificial compressibilitiy
    %keyboard;
    lhs3=p+dt*(c02*Div2D(u,v)-LIFT*(Fscale.*(c02*nx.*du+c02*ny.*dv -lambda.*dp)/2));
    %lhs3= (c02*Div2D(u,v)-LIFT*(Fscale.*(c02*nx.*du+c02*ny.*dv -lambda.*dp)/2));
    %lhs3=Div2D(u,v);
    
    result=[lhs1(:);lhs2(:);lhs3(:)];
     
end