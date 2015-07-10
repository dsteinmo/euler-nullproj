function [Q] = Euler2D(Q, FinalTime, BC)

% function [Q] = Euler2D(Q,FinalTime,BC);
% Purpose  : Integrate 2D Euler equations using a 3rd order SSP-RK

Globals2D;

% Initialize filter
Filt = CutOffFilter2D(N,0.95);

% compute initial timestep
gamma = 1.4;
dt = EulerDT2D(Q, gamma);
time = 0; tstep=1;

% storage for low storage RK time stepping
rhsQ = 0*Q; resQ = 0*Q;

% filter initial solution
for n=1:4, Q(:,:,n) = Filt*Q(:,:,n); end;

%Test 2D quadrature integration
%myfcn = sin(2*pi*x/10).^2 + 1;
% myfcn = (sin(2*pi*x/10).^2).*(cos(2*pi*y/10).^2);
% 
% myfcnmodes = V\myfcn;  %transform nodes into modes
% modezero = myfcnmodes(1,:);  %Get mode zero on each element
% Po = 1/sqrt(2); %Scaling factor, since the height of order 0 ortho. legendre poly is 1/sqrt(2)
% Kscale = 2*J(1,:);
% integral = Po*(modezero*Kscale') %dot product to sum up integral on each element


figure(3);
%%PlotField2D(N,x,y,myfcn); view([0 90]); axis([0 10 -5 5]); drawnow; colorbar;
%keyboard;

% outer time step loop 
figure(1);
while (time<FinalTime)

  % check to see if we need to adjust for final time step
  if(time+dt>FinalTime)
    dt = FinalTime-time;
  end
  
  for INTRK = 1:5    
    % compute right hand side of compressible Euler equations
    rhsQ  = EulerRHS2D(Q, time, BC);
    
    % filter residual
    for n=1:4, rhsQ(:,:,n) = Filt*rhsQ(:,:,n); end;
    
    % initiate and increment Runge-Kutta residuals
    resQ = rk4a(INTRK)*resQ + dt*rhsQ;  
    
    % update fields
    Q = Q+rk4b(INTRK)*resQ;  
  end;
  
  % Increment time and compute new timestep
  time = time+dt;
  dt = EulerDT2D(Q, gamma);
  
  if mod(tstep,100)==0 || time >= FinalTime
    PlotField2D(N, x, y, Q(:,:,4));
    view([-86 38])
    title(['t=' num2str(time)]); 
    xlabel('x');
    ylabel('y');
    zlabel('z');
    colorbar; drawnow;
  end
  tstep = tstep+1;
end;
return;
