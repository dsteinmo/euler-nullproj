

%takes in predicted (u,v) as input, returns projected (u,v,p)
%assumes bordering trick
function result = pressproj_hnd(RHS,dt,ll,uu,pp,qq)
    %Globals2D;
    global Np K Nfaces Nfp vmapM vmapP LIFT Fscale nx ny mapW vmapW MassMatrix J
    global EToE EToF BCType Wall N x y
    
    postprocessop=false;
    
%     size(RHS)
%     keyboard;
    ustar=RHS(1:Np*K); ustar=reshape(ustar,Np,K);
    vstar=RHS(Np*K+1:2*Np*K); vstar=reshape(vstar,Np,K);
    
    
    %compute divergence of u*
  divu = Div2D(ustar,vstar);  %compute weak divergence
  du = zeros(Nfaces*Nfp,K); dv = zeros(Nfaces*Nfp,K);
  du(:) = (ustar(vmapM)-ustar(vmapP));  %form field differences at interfaces
  dv(:) = (vstar(vmapM)-vstar(vmapP));
  fluxu = (nx.*du + ny.*dv)/2.0; %compute central flux
  preRHS = (1/dt)*(divu - LIFT*(Fscale.*fluxu)); %'full' divergence
  
  RHS = MassMatrix*(J.*preRHS);

  %Evaluate BC contribution to RHS
  qW = zeros(Nfp*Nfaces, K);
  %Neumann BC = (1/dt)*(n . ustar)
  qW(mapW) = (1/dt)*(nx(mapW).*ustar(vmapW) + ny(mapW).*vstar(vmapW)); %this is necessary and right.
  Aqbc = PoissonIPDGbc2Dquad([],qW);
  
  RHS = RHS(:) - Aqbc(:); %think this is right
  
  p=qq*(uu\(ll\(pp*[RHS;0])));
  p = p(1:end-1);
  
  p = reshape(p,Np,K);
  
  [px,py] = Grad2D(p);

  %time-step pressure
  dp = zeros(Nfaces*Nfp,K);
  dp(:) = p(vmapM)-p(vmapP);
  u = ustar - dt*(px - LIFT*(Fscale.*(nx.*dp/2)));   %choice of flux
  v = vstar - dt*(py - LIFT*(Fscale.*(ny.*dp/2)));   %doesn't matter much here (can go fully interior)
  %u = ustar - dt*px;
  %v = vstar - dt*py;
  
%   uraw=u;
%   vraw=v;
  
  if postprocessop==true
      disp('hi');
      for k1=1:K
          for f1=1:Nfaces

              k2 = EToE(k1,f1); f2 = EToF(k1,f1);
              fidM  = (k1-1)*Nfp*Nfaces + (f1-1)*Nfp + (1:Nfp);
              vidM = vmapM(fidM); Fm1 = mod(vidM-1,Np)+1;
              vidP = vmapP(fidM); Fm2 = mod(vidP-1,Np)+1;

              id = 1+(f1-1)*Nfp + (k1-1)*Nfp*Nfaces;
              lnx = nx(id);  lny = ny(id); %lsJ = sJ(id);


              %keyboard;
              u1 = u(Fm1,k1); u2=u(Fm2,k2);
              v1 = v(Fm1,k1); v2=v(Fm2,k2);

              %ustar1 = ustar(Fm1,k1); %get velocity along this face
              %vstar1 = vstar(Fm1,k1);
              %ustar2 = ustar(Fm2,k2);
              %figure(1);
              %plot(1:Nfp,ustar1,1:Nfp,ustar2);
              %keyboard;

              %note that this assumes (lnx,lny) are scalars

              %ustar1dotn = mean(lnx*ustar1+lny*vstar1); %trying this.

              un1 = lnx*u1 + lny*v1; 
              ut1 =-lny*u1 + lnx*v1; 
              un2 = lnx*u2 + lny*v2; 
              ut2 =-lny*u2 + lnx*v2; 

              %pick appropriate side to glue from based on upwinding
              if BCType(k1,f1) == Wall
                  %strongly impose no normal flow
                  un=0*u1;
              else
                  %un = (ustar1dotn>=0)*un1 + (ustar1dotn<0)*un2; %this is right
                  %un = (ustar1dotn<0)*un1 + (ustar1dotn>=0)*un2; %this is wrong
                  un = 0.5*(un1+un2);
              end

              un1 = un;
              %un2 = un; %seeing if I need this part.

              %rotate back
              u1 = lnx.*un1 - lny.*ut1;
              v1 = lny.*un1 + lnx.*ut1;
              
              %u2 = lnx.*un2 - lny.*ut2;
              %v2 = lny.*un2 + lnx.*ut2;

              u(Fm1,k1)=u1; 
              v(Fm1,k1)=v1; 
              
              %u(Fm2,k2)=u2;
              %v(Fm2,k2)=v2;

                %like spectral elements -- bad
                %u(Fm1,k1)=.5*(u1+u2);
                %v(Fm1,k1)=.5*(v1+v2);
                
                %u(Fm2,k2)=.5*(u1+u2);
                %v(Fm2,k2)=.5*(v1+v2);

          end
      end
  end
  
%  figure(1);
%  subplot(2,2,1);
%  pf2d(N,x,y,ustar); shading interp; colorbar;
%  subplot(2,2,2);
%  pf2d(N,x,y,uraw); shading interp; colorbar;
%  subplot(2,2,3);
%  pf2d(N,x,y,u); shading interp; colorbar;
%  subplot(2,2,4);
%  pf2d(N,x,y,u-uraw); shading interp; colorbar;
%  drawnow;
%  
%  figure(2);
%   subplot(2,2,1);
%  pf2d(N,x,y,vstar); shading interp; colorbar;
%  subplot(2,2,2);
%  pf2d(N,x,y,vraw); shading interp; colorbar;
%  subplot(2,2,3);
%  pf2d(N,x,y,v); shading interp; colorbar;
%  subplot(2,2,4);
%  pf2d(N,x,y,v-vraw); shading interp; colorbar;
%  
%  drawnow;
%  
%  keyboard;
  
  result=[u(:);v(:);p(:)];
     
end