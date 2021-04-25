function yi = myinterp1_equi(x,y,xi)

% Usage: yi = myinterp1_equi(x,y,xi)
% This function is very fast also because it does not do any input check:
% it assumes that you know what you are doing!
% INPUTS:
% x and y are N*1 vectors
% x MUST be monotonically increasing (grid) AND equispaced
% xi is the query point and it MUST be a scalar
% n is the length of vector x and y
% OUTPUT
% yi is the value of the interpolated function at query point xi 

% -------------------------------------------------------------------------
% Copyright © 2018 by Alessandro Di Nola. All rights 
% reserved.
% -------------------------------------------------------------------------

n = size(x,1);

% xi is between x(j) and xx(j+1)
% j = 0 and j = n means x is out of range!!!!!!!!!!!!!!
% equispaced grid
step = x(2)-x(1);
j = max(min(ceil((xi-x(1))/step),n-1),1);

slope = (y(j+1)-y(j))/(x(j+1)-x(j));
yi = y(j)+(xi-x(j))*slope;

end