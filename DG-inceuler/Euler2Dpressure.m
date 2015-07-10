
function [Q] = Euler2Dpressure(Qn, FinalTime, BC,g)

% function [Q] = Euler2D(Q,FinalTime,BC);
% Purpose  : Integrate 2D Euler equations using a 3rd order SSP-RK

Globals2D;

Lap = -PoissonIPDG2Dquad; 
OP = [Lap ones(Np*K,1);
     ones(1,Np*K) 0];   %linear algebra trick to get rid of mean (i.e., the null eigenmode)

[ll,uu,pp,qq] = lu(OP);

clear OP

% Initialize filter - be careful not to 'over filter'
Filt = Filter2Dquad(N,0.9*N,4);%order 4 runs

FILTER_TRACER_RHS = false; %seem to get away with 'false'

H = abs(max(y(:))-min(y(:)));

% compute initial timestep -- time-stepping is not adaptive (yet).
CFL=0.85;
dt=CFL*EulerDT2D(Qn,g,H);

tstep=1;
time(tstep)=0;
maxdivseries(tstep)=0;

% storage for AB3 time stepping
rhsQnm1 = 0*Qn;
rhsQnm2 = 0*Qn;

Qn(:,:,4) = 0*Qn(:,:,4);  %should get rid of 4th index, it's not used.
%dt=1e-2;
Qnp1 = zeros(size(Qn));
bn = ab3coefs([time(1)+dt time(1) time(1) time(1)]);

%build data necessary for element-wise projection
projectionMatrix = buildProjectionMatrixquad(N,r,s,V);
rotData = buildCoordRotationData_quad;

%Project initial velocity (NEW!) to ensure div-free.
u=Qn(:,:,2); v=Qn(:,:,3);
[u,v] = globalDivFreeProjection(u,v,projectionMatrix,rotData);
Qn(:,:,2)=u;
Qn(:,:,3)=v;

%output initial data
Qnp1=Qn;
save(sprintf('out%07d.mat',0),'Qnp1','N','K','x','y','time');

disp('entering main time loop...');
figure(1); clf;
while (time(tstep)<FinalTime)
  
   %get advective RHS's
   rhsQn = EulerRHS2D(Qn,time,BC,g);  
  
   if FILTER_TRACER_RHS == true
       rhsQn(:,:,1) = Filt*rhsQn(:,:,1);
   end
   
   %advective step (AB3)
   Qnp1 = Qn + bn(1)*rhsQn + bn(2)*rhsQnm1 + bn(3)*rhsQnm2;
 
   %filter velocity fields after advective step (a must!)
   for n=2:3
       Qnp1(:,:,n) = Filt*Qnp1(:,:,n);
   end
   
   %get predicted velocity
   ustar = Qnp1(:,:,2); vstar = Qnp1(:,:,3); %for reg projection
   
   RHS = [ustar(:); vstar(:)];
   result = pressproj_hnd(RHS,dt,ll,uu,pp,qq);
   u = reshape(result(1:end/3),Np,K);
   v = reshape(result(end/3+1:2*end/3),Np,K);
   p = reshape(result(2*end/3+1:end),Np,K);
   
   [px,py] = Grad2D(p); %for diagnostic only
   
   %now do post-processing projection to get a truly
   %div free field
   [u,v] = globalDivFreeProjection(u,v,projectionMatrix,rotData);
   
   %put corrected velocities back into 'packed' form.
   Qnp1(:,:,2)=u;
   Qnp1(:,:,3)=v;
   
   % Increment time and compute new timestep
   tstep = tstep+1;
   time(tstep)=time(tstep-1)+dt;
  
  %compute new time-step
  if tstep == 2
      %AB2
      bn=ab3coefs([time(tstep)+dt,time(tstep),time(tstep-1),time(tstep-1)]);
  else
      bn=ab3coefs([time(tstep)+dt,time(tstep),time(tstep-1),time(tstep-2)]);
  end
  
  divu = Div2D(u,v);
  maxdivseries(tstep)=max(abs(divu(:)));
  divnormseries(tstep) = sqrt(dgintquad(divu.^2,V,J));
  
 if mod(tstep,100)==0 || time(tstep) >= FinalTime ||  tstep==2
    

    
    clf;
    subplot(2,3,1);
    pf2dquad(N, x, y, Qnp1(:,:,1)+1); %go back to "full"
    title(['\rho @ t=' num2str(time(tstep))]); 
    xlabel('x');
    ylabel('z');
    %caxis([1 1.003]);
    colorbar
    
    
    subplot(2,3,2);
    pf2dquad(N,x,y,px);
    title('px');
    axis tight;
    colorbar;
    
    subplot(2,3,3);
    pf2dquad(N,x,y,py);
    title('py');
    axis tight;
    colorbar;
    
    subplot(2,3,4);
    %pf2d(N,x,y,LIFT*((u(vmapM)-u(vmapP))./(u(vmapM)+u(vmapP)))); colorbar;
    pf2dquad(N,x,y,Qnp1(:,:,2)); colorbar;
    subplot(2,3,5);
    %pf2d(N,x,y,LIFT*((v(vmapM)-v(vmapP))./(v(vmapM)+v(vmapP)))); colorbar;
    pf2dquad(N,x,y,Qnp1(:,:,3)); colorbar;
    %subplot(2,3,6);
    %pf2d(N,x,y,divu); colorbar;
    %title('div u');
    subplot(2,3,6);
    p_hyd = -g*Qnp1(:,:,1).*y;
    p_hyd = p_hyd - mean(p_hyd(:));
    
    pf2dquad(N,x,y,p-p_hyd);
    title('pressure perturbation'); colorbar;
    drawnow;
    %bn
    
    save(sprintf('out%07d.mat',tstep),'x','y','Np','K','N','Qnp1','p','time','tstep','maxdivseries','divu','divnormseries');
    disp(['outputting at t=' num2str(time(tstep))]);
    
    
 end %end output loop
 % keyboard;
  
 if any(isnan(Qnp1(:)))
     disp('NaN''s!');
     Q=Qnp1;
     return
 end


  %rotate fields for next time-step
  Qn=Qnp1;
  rhsQnm2=rhsQnm1;
  rhsQnm1=rhsQn;
end;
Q=Qnp1;
return;
