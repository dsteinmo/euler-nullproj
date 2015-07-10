function [IM] = InterpMatrix2Dquad(rout, sout)

% function [IM] = InterpMatrix2D(rout, sout)
% Purpose: Compute local elemental interpolation matrix

Globals2D;
 
% compute Vandermonde at (rout,sout)
Vout = Vandermonde2Dquad(N, rout, sout);

% build interpolation matrix
IM = Vout*invV;
return

  
  
