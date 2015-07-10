%unit test for bandLimtFourierInterp2D - modified from:
% Trefethen, "Spectral Methods in MATLAB", p3.m.
clear
close all;

%make sure that the functions we're calling are visible to matlab
addpath ../

%do you want to plot results?
GRAPHICAL_OUTPUT = true;

%set error tolerance
errtol = 1e-4;

%build grid periodic
xmax = 10; 
xmin = -xmax;
Lx=(xmax-xmin);
Nx=40;
h = Lx/Nx;
x=(1:Nx)*h-xmax;

%pick some interpolation points
xout = randn(1,200)*xmax;
xout = xout(xout <= x(end) & xout >= x(1)); %make sure we're _interpolating_

%test #1
%set signal to interpolate, and do interpolation to above points
v = sech(x).^2;
p = bandLimFourierInterp1D(x,v,xout);

%check error - should get better if we refine the grid
err = norm(sech(xout(:)).^2-p(:),2)/length(p);
if err < errtol
    disp(['Test of bandLimFourierInterp1D.m (decaying signal) PASSED with err=' num2str(err)]);
else
    disp(['Test of bandLimFourierInterp1D.m (decaying signal) FAILED with err=' num2str(err)]);
end

%%plot results if you want
if GRAPHICAL_OUTPUT == true
    figure(1); clf;
    plot(x,v,'.-',xout,p,'.r'); 
    legend('original signal','interpolated values');
    axis([xmin xmax -.1 1.1]); 
end

%test #2
%set signal to interpolate, and do interpolation to below points
xout = xmin:(xmax-xmin)/200:xmax;
xout = xout(xout <= x(end) & xout >= x(1)); %make sure we're _interpolating_
v = sin(2*pi*x/xmax);
p = bandLimFourierInterp1D(x,v,xout);

err = norm(sin(2*pi*xout(:)/xmax)-p(:),2)/length(p);
if err < errtol
    disp(['Test of bandLimFourierInterp1D.m (odd signal) PASSED with err=' num2str(err)]);
else
    disp(['Test of bandLimFourierInterp1D.m (odd signal) FAILED with err=' num2str(err)]);
end

%%plot results if you want
if GRAPHICAL_OUTPUT == true
    figure(2); clf;
    plot(x,v,'.-',xout,p,'.r'); 
    legend('original signal','interpolated values','location','southwest');
end

%test #3
%set signal to interpolate, and do interpolation to below points
xout = xmin:(xmax-xmin)/200:xmax;
xout = xout(xout <= x(end) & xout >= x(1)); %make sure we're _interpolating_
v = cos(2*pi*x/xmax);
p = bandLimFourierInterp1D(x,v,xout);

err = norm(cos(2*pi*xout(:)/xmax)-p(:),2)/length(p);
if err < errtol
    disp(['Test of bandLimFourierInterp1D.m (even signal) PASSED with err=' num2str(err)]);
else
    disp(['Test of bandLimFourierInterp1D.m (even signal) FAILED with err=' num2str(err)]);
end

%%plot results if you want
if GRAPHICAL_OUTPUT == true
    figure(3); clf;
    plot(x,v,'.-',xout,p,'.r'); 
    legend('original signal','interpolated values','location','southwest');
end

%%test 4 - an aperiodic signal -> this is expected to fail due to gibbs
%%oscillations, since input data must be periodic
% v = tanh(x);
% xout=x+h/2;
% xout = xout(xout <= x(end) & xout >= x(1)); %make sure we're _interpolating_
% p = bandLimFourierInterp1D(x,v,xout);
% 
% err = norm(tanh(xout(:))-p(:),2)/length(p);
% if err < errtol
%     disp(['Test of bandLimFourierInterp1D.m (aperiodic) PASSED (somehow?!) with err=' num2str(err)]);
% else
%     disp(['Test of bandLimFourierInterp1D.m (aperiodic) FAILED (as expected) with err=' num2str(err)]);
% end
% 
% figure(4); clf;
% plot(x,v,'.-',xout,p,'.r'); 
% legend('original signal','interpolated values');
% axis([xmin xmax -1.1 1.1]); 