
function [Q] = Euler2Dextrap(Qn, FinalTime, BC,g)

% function [Q] = Euler2D(Q,FinalTime,BC);
% Purpose  : Integrate 2D Euler equations using a 3rd order SSP-RK

Globals2D;

%Build poisson operator
Lap = -PoissonIPDG2Dquad; 
OP = [Lap ones(Np*K,1);
     ones(1,Np*K) 0];   %linear algebra trick to get rid of mean (i.e., the null eigenmode)

[Pll,Puu,Ppp,Pqq] = lu(OP);

% Build misfit operators.
C0 = build_C0_misfit();
B0 = build_udotn_misfit();
C1x = build_c1x_misfit();
C1y = build_c1y_misfit();

% Build C0 continuity misfit for velocity field.
C0_vel = [C0+C1y C1x; C1y C0+C1x];

% Combine continuity and boundary conditions into single operator for velocity field.
E = C0_vel + B0;

% Build basis for the nullspace of the divergence operator.
Nd = build_div_nullspace();

% Form normal equation.
tau = 0;
T = (Nd'*Nd) + tau * ( E * Nd )' * E * Nd;

%Pre-factorize normal operator
[ll,uu,pp,qq] = lu(T);

% Initialize filter - be careful not to 'over filter'
Filt = Filter2Dquad(N,0.9*N,4);%order 4 runs

FILTER_TRACER_RHS = false; %seem to get away with 'false'

H = abs(max(y(:))-min(y(:)));

% compute initial timestep -- time-stepping is not adaptive (yet).
CFL=0.9;
dt=CFL*EulerDT2D(Qn,g,H);

tstep=1;
time(tstep)=0;
maxdivseries(tstep)=0;

% storage for AB3 time stepping
rhsQnm1 = 0*Qn;
rhsQnm2 = 0*Qn;

Qn(:,:,4) = 0*Qn(:,:,4);  %should get rid of 4th index, it's not used.

Qnp1 = zeros(size(Qn));
%build data necessary for element-wise projection
projectionMatrix = buildProjectionMatrixquad(N,r,s,V);
rotData = buildCoordRotationData_quad;

%Project initial velocity (NEW!) to ensure div-free.
%u=Qn(:,:,2); v=Qn(:,:,3);
%[u,v] = globalDivFreeProjection(u,v,projectionMatrix,rotData);
%Qn(:,:,2)=u;
%Qn(:,:,3)=v;

%output initial data
Qnp1=Qn;
save(sprintf('out%07d.mat',0),'Qnp1','N','K','x','y','time');


px_nmh = zeros(Np,K);
py_nmh = zeros(Np,K);
gradp_nmh = zeros(Np,K,4);

disp('entering main time loop...');
figure(1); clf;
while (time(tstep)<FinalTime)
  
   %get advective RHS's
   rhsQn = EulerRHS2D(Qn,time,BC,g);  
  
   if FILTER_TRACER_RHS == true
       rhsQn(:,:,1) = Filt*rhsQn(:,:,1);
   end
   
   %advective step (predictor)
   Qnp1_estimate = Qn + dt*rhsQn - dt*gradp_nmh;

   
   %Form estimate at n+1/2 step.
   Qnph = 0.5*(Qnp1_estimate + Qn);

   %Get estimate for convective/source terms at half step.
   rhsQnph = EulerRHS2D(Qnph,time,BC,g);

   %extrapolate to get corrected estimate for grad p at n+1/2 step.
   gradp_nph = rhsQnph - rhsQn + gradp_nmh;

   %advective step (corrector)
   Qnp1 = Qn + dt*rhsQnph - dt*gradp_nph;


   %filter velocity fields after advective step (a must!)
   for n=2:3
       Qnp1(:,:,n) = Filt*Qnp1(:,:,n);
   end
   
   %get predicted velocity
   ustar = Qnp1(:,:,2); vstar = Qnp1(:,:,3); %for reg projection

   RHS = [ustar(:); vstar(:)];
   result = pressproj_hnd(RHS,dt,Pll,Puu,Ppp,Pqq);
   u = reshape(result(1:end/3),Np,K);
   v = reshape(result(end/3+1:2*end/3),Np,K);
   p = reshape(result(2*end/3+1:end),Np,K);

   enerpre = 0.5*(u.^2 + v.^2);

   uvec = [u(:); v(:)];
  
   % Set up the right-hand-side for the normal equation.
   RHS = Nd' * [ uvec ];

   % Solve for new velocity with LUPQ factors
   lambda=qq*(uu\(ll\(pp*RHS)));

   % Construct the updated velocity.
   unew = Nd * lambda;

   %DEREK: remove comments before pushing. - with this in, it died at t=145. without this, dies at t=11. without this and filtering: t=200 or more?
   u = reshape(unew(1:end/2),Np,K);
   v = reshape(unew(end/2+1:end),Np,K);


   %u = Filt*u;
   %v = Filt*v;
   enerpost = 0.5*(u.^2 + v.^2);

   Qnp1(:,:,2) = u;
   Qnp1(:,:,3) = v; 

   gradp_nph = (1/dt)*(Qnp1 - Qn);
   gradp_nph(:,:,1) = zeros(Np,K);
   px = gradp_nph(:,:,2);
   py = gradp_nph(:,:,3);

   %compute divergence (for diagnostic).
   divu = Div2D(u,v);

   %put corrected velocities back into 'packed' form.
   Qnp1(:,:,2)=u;
   Qnp1(:,:,3)=v;
   
   % Increment time and compute new timestep
   tstep = tstep+1;
   time(tstep)=time(tstep-1)+dt;
  
  
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
    title('pz');
    axis tight;
    colorbar;
    
    subplot(2,3,4);
    %pf2d(N,x,y,LIFT*((u(vmapM)-u(vmapP))./(u(vmapM)+u(vmapP)))); colorbar;
    pf2dquad(N,x,y,Qnp1(:,:,2)); colorbar;
    title('u');
    subplot(2,3,5);
    %pf2d(N,x,y,LIFT*((v(vmapM)-v(vmapP))./(v(vmapM)+v(vmapP)))); colorbar;
    pf2dquad(N,x,y,Qnp1(:,:,3)); colorbar;
    title('w');
    subplot(2,3,6);
    pf2dquad(N,x,y,divu); colorbar;
    title('div u');
    %subplot(2,3,6);
    %p_hyd = -g*Qnp1(:,:,1).*y; 
    %pf2dquad(N,x,y,enerpost-enerpre);
    %title('energy difference'); colorbar;
    drawnow;
    %bn
    
    save(sprintf('out%07d.mat',tstep),'x','y','Np','K','N','Qnp1','time','tstep','maxdivseries','divu','divnormseries');
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
  gradp_nmh = gradp_nph;
  %rhsQnm2=rhsQnm1;
  %rhsQnm1=rhsQn;
end;
Q=Qnp1;
return;
