%unit test for bandLimtFourierInterp2D
clear;
close all;

%make sure that the functions we're calling are visible to matlab
addpath ../ 

%do you want to plot results?
GRAPHICAL_OUTPUT = true;

%set error tolerance
errtol = 1e-3;

%set optional max memory argument
maxMem=[];
%maxMem=4e8;

%build structured 64x64 grid
Nx=64;
Ny=64;
dx = 2*pi/Nx;
dy = 2*pi/Ny;
[x,y]= meshgrid((1:Nx)*dx,(1:Ny)*dy);
f=sin(x).*sin(y); %function to sample

%construct set of points to interpolate at.
xout= randn(100,1)*pi + pi;
yout= randn(100,1)*pi + pi;
xoutnew = xout(xout >x(1,1) & xout < x(end,end) & yout >y(1,1) & yout < y(end,end));
yout = yout(xout >x(1,1) & xout < x(end,end) & yout >y(1,1) & yout < y(end,end));
xout =xoutnew; clear xoutnew;

fout = bandLimFourierInterp2D(x,y,f,xout,yout,maxMem);

%check error - should get better if we refine the grid
err = norm(sin(xout(:)).*sin(yout(:))-fout(:),2)/length(xout);
if err < errtol
    disp(['Test of bandLimFourierInterp2D.m PASSED with err=' num2str(err)]);
else
    disp(['Test of bandLimFourierInterp2D.m FAILED with err=' num2str(err)]);
end

%%plot results, if you want
if GRAPHICAL_OUTPUT == true
    figure(1); clf;
    surf(x,y,f); shading interp; colorbar; drawnow;
    hold on;
    plot3(xout,yout,fout,'*r')
end

%test 2: test tensor-product output functionality.
[xxout,yyout] = meshgrid(sort(xout,'ascend'),sort(yout,'ascend'));
fout_tensor = bandLimFourierInterp2D(x,y,f,xxout,yyout,maxMem);

%did we reshape properly?
if all(size(fout_tensor) == size(xxout))
    disp('Test of bandLimFourierInterp2D.m (tensor product output) PASSED');
else
    disp('Test of bandLimFourierInterp2D.m (tensor product output) FAILED');
end

%%plot results, if you want
if GRAPHICAL_OUTPUT == true
    figure(2); clf;
    surf(xxout,yyout,fout_tensor); shading interp; colorbar; drawnow;
    hold on;
    plot3(xxout,yyout,fout_tensor,'*r');
end

%test 3: test that we're getting boundaries right (odd/odd symmetry)
Nout=64;
[xxout,yyout] = meshgrid(linspace(dx,2*pi,Nout),linspace(dy,2*pi,Nout));
fout_tensor = bandLimFourierInterp2D(x,y,f,xxout,yyout,maxMem);

%check error - should get better if we refine the grid
err = norm(sin(xxout(:)).*sin(yyout(:))-fout_tensor(:),2)/length(xout);
if err < errtol
    disp(['Test of bandLimFourierInterp2D.m (odd/odd boundary test) PASSED with err=' num2str(err)]);
else
    disp(['Test of bandLimFourierInterp2D.m (odd/odd boundary test) FAILED with err=' num2str(err)]);
end

%%plot results, if you want
if GRAPHICAL_OUTPUT == true
    figure(3); clf;
    surf(xxout,yyout,fout_tensor); shading interp; colorbar; drawnow;
    hold on;
    plot3(xxout,yyout,fout_tensor,'*r');
end

%test 4: test that we're getting boundaries right (odd/even symmetry)
f=sin(x).*cos(y); %function to sample
Nout=64;
[xxout,yyout] = meshgrid(linspace(dx,2*pi,Nout),linspace(dy,2*pi,Nout));
fout_tensor = bandLimFourierInterp2D(x,y,f,xxout,yyout,maxMem);

%check error - should get better if we refine the grid
err = norm(sin(xxout(:)).*cos(yyout(:))-fout_tensor(:),2)/length(xout);
if err < errtol
    disp(['Test of bandLimFourierInterp2D.m (odd/even boundary test) PASSED with err=' num2str(err)]);
else
    disp(['Test of bandLimFourierInterp2D.m (odd/even boundary test) FAILED with err=' num2str(err)]);
end

%%plot results, if you want
if GRAPHICAL_OUTPUT == true
    figure(4); clf;
    surf(xxout,yyout,fout_tensor); shading interp; colorbar; drawnow;
    hold on;
    plot3(xxout,yyout,fout_tensor,'*r');
end

%test 5: test that we're getting boundaries right (odd/even symmetry)
f=cos(x).*cos(y); %function to sample
Nout=64;
[xxout,yyout] = meshgrid(linspace(dx,2*pi,Nout),linspace(dy,2*pi,Nout));
fout_tensor = bandLimFourierInterp2D(x,y,f,xxout,yyout,maxMem);

%check error - should get better if we refine the grid
err = norm(cos(xxout(:)).*cos(yyout(:))-fout_tensor(:),2)/length(xout);
if err < errtol
    disp(['Test of bandLimFourierInterp2D.m (even/even boundary test) PASSED with err=' num2str(err)]);
else
    disp(['Test of bandLimFourierInterp2D.m (even/even boundary test) FAILED with err=' num2str(err)]);
end

%%plot results, if you want
if GRAPHICAL_OUTPUT == true
    figure(5); clf;
    surf(xxout,yyout,fout_tensor); shading interp; colorbar; drawnow;
    hold on;
    plot3(xxout,yyout,fout_tensor,'*r');
end