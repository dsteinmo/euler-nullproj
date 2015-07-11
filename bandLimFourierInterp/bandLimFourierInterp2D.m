%BANDLIMFOURIERINTERP2D 2D Band-Limited Fourier Interpolation
%   fout = bandLimFourierInterp2D(x,y,f,xout,yout) takes periodic signal 
%   f sampled at equispaced grid (x,y) and interpolates to arbitrary set of
%   points (xout,yout). [x,y] should be a tensor-product grid generated with
%   meshgrid. [xout,yout] can either be a tensor-product grid or a list
%   of points where xout and yout are each 1D arrays.
%
%   fout = bandLimFourierInterp2D(x,y,f,xout,yout,maxMem) is the same as
%   above but with optional argument maxMem that allows the specification of 
%   maximum amount of memory (in bytes) allowed to be allocated by temporary 
%   3D tensor-product arrays. Default is 200 MB.
%
%   Uses band-limited interpolation formula in tensor product form to avoid
%   'for' loops. cf. Trefethen, "Spectral Methods in MATLAB", p3.m.
%
%   Should run in both MATLAB and GNU Octave, please report any bugs to 
%   the author: dsteinmo "at" uwaterloo.ca
%
%   Author: Derek Steinmoeller, University of Waterloo. Copyright (C) 2012.
%   Improvements by: Mike Dunphy (improved speed/memory efficiency)
%                    Chris Subich (help with correcting periodicity bugs)
%
%   This program is free software; you can redistribute it and/or
%   modify it under the terms of the GNU General Public License version 3
%   as published by the Free Software Foundation.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA
function fout =  bandLimFourierInterp2D(x,y,f,xout,yout,maxMem)

    %deal with optional input argument
    MAXMEMDEFAULT = 209715200; %default = 200 MB
    
    if ~exist('maxMem','var')
        maxMem = MAXMEMDEFAULT;
    else
        if isempty(maxMem)
            maxMem = MAXMEMDEFAULT;
        end
        %else, input is probably acceptable
    end
    
    %get 1D versions of x & y-grids
    x1d = x(1,:);
    y1d = y(:,1);
    Nx = length(x1d);
    Ny = length(y1d);

    %get grid-spacing
    dx = abs(x1d(2)-x1d(1));
    dy = abs(y1d(2)-y1d(1));
    
    %get domain length & width
    Lx = dx*length(x);
    Ly = dy*length(y);
    
    %depending on the parity of the original signal's number of nodes (in
    %each direction) choose the appropriate periodic sinc.
    if mod(Nx,2)==1
        psincX = @(xx) bsinc(pi*xx,dx,Lx);
    else
        psincX = @(xx) bsinc(pi*xx,dx,Lx).*cos(pi*xx/Lx);
    end
    
    if mod(Ny,2)==1
        psincY = @(yy) bsinc(pi*yy,dy,Ly);
    else
        psincY = @(yy) bsinc(pi*yy,dy,Ly).*cos(pi*yy/Ly);
    end

    %determine if we need to return a 1D or 2D array.
    szout = size(xout);
    if szout(1) > 1 && szout(2)> 1
        TENSOR_PROD_OUTPUT = true;
    else
        TENSOR_PROD_OUTPUT = false;
    end

    %make output points into column vectors for rest of calculations
    xout=xout(:);
    yout=yout(:);
    Nout=length(xout);

    %calculate acceptable size of 3D arrays
    sizeOf3D =maxMem/2;               %Calculation needs 2 3D arrays
    maxEntries3D = sizeOf3D/8;        %populated by 8-byte entries.
    %so third dimension should be at most this big:
    zlen = floor(maxEntries3D/Nx/Ny);

    if zlen == 0
        error('Input data too big or maxMem too small. Try making maxMem larger (default is 209715200 bytes).');
    end

    %if we can fit the entire calculation in memory
    %at the same time, then proceed as normal (as in older version)
    if Nout < zlen
        fout = doSincInterpFast(f,psincX,psincY,x1d,y1d,xout,yout);
    else
        %if we can't fit it into memory, split it up into blocks of length
        %zlen and process them sequentially with a for-loop.
        Nblocks = ceil(Nout/zlen);
        fout = zeros(Nout,1);

        for n=1:Nblocks
            if n==Nblocks
                %last block could be smaller
                range = (zlen*(n-1)+1):Nout;
            else
                range = (zlen*(n-1)+1):zlen*n;
            end

            fout(range) = doSincInterpFast(f,psincX,psincY,x1d,y1d,xout(range),yout(range));
        end
    end

    %reshape data, if we're supposed to
    if TENSOR_PROD_OUTPUT == true
        fout = reshape(fout,szout(1),szout(2));
    end
end

%Helper function that actually does the interpolation - to be called repeatedly
function fout = doSincInterpFast(f,sincX,sincY,x1d,y1d,xout,yout)
    Nx = length(x1d); Ny = length(y1d); Nout = length(xout);

    fprintf('Working on %d points...\n',Nout);

    %evaluate all the shifted sinc functions using bsxfun
    x1d=reshape(x1d,[1 Nx 1]);      % transpose
    xout=reshape(xout,[1 1 Nout]);  % appropriately.
    shiftX = bsxfun(@minus,xout,x1d);  % result is 2D array

    y1d=reshape(y1d,[Ny 1 1]);      % transpose
    yout=reshape(yout,[1 1 Nout]);  % appropriately.
    shiftY = bsxfun(@minus,yout,y1d);  % result is 2D array
    
    %evaluate kernel function and calculate the contribution
    %due to each grid point
    %kernel = bsxfun(@times,sincX(shiftX),sincY(shiftY)); %3D array
    %fout = bsxfun(@times,f,kernel); %3D array
    %...this step can be optimized by combining into one line:
    fout = bsxfun(@times,f,bsxfun(@times,sincX(shiftX),sincY(shiftY)));

    %double-contract (sum over) them all to 
    %get interpolated data at (xout,yout).
    fout = squeeze(sum(sum(fout,2),1));
end
% Helper function: Baseline sinc -- returns sin(x/dx)/sin(x/L)*(dx/L), since
% the special case x == 0 gives a naive implementation NaN issue
function z = bsinc(x,dx,L)
    z = zeros(size(x));
    z(x==0) = 1;
    z(x~=0) = sin(x(x~=0)/dx)./sin(x(x~=0)/L)*(dx/L);
    return
end