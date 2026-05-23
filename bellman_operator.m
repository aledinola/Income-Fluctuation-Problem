function [Tv,c_pol] = bellman_operator(v,par)
% This function applies the BO to the function v, returning Tv
% I tested several interpolation methods (matlab built-in, my interpolation
% function myinterp1_v1, etc.
% Author: Alessandro Di Nola

%Unpack inputs
R          = par.R;
PZ         = par.PZ;
beta       = par.beta;
b          = par.b;
asset_grid = par.asset_grid;
z_vals     = par.z_vals;
na         = par.na;
nz         = par.nz;
tiny       = par.tiny;
u          = par.u;
tol_gss    = par.tol_gss;

% Total number of points for state variables
nx = na*nz;

%Solve for RHS of Bellman equation
%{
TV(a,z) = max{ u(c)+beta*E[V(R*a+z-c,z')|z] }
where the max is over c in the set [0,R*a+z+b]
Clearly, a'=R*a+z-c
%}

c_pol = zeros(nx,1);
Tv    = zeros(nx,1);

parfor ix=1:nx
    [ia,iz] = ind2sub([na,nz],ix);
    z_today = z_vals(iz);

    % Compute continuation value
    PZ_iz = PZ(iz,:);
    %Now sum over all z'
    EVz = v*PZ_iz'; %this is EV(a'), given z

    % Do the maximization
    a_today  = asset_grid(ia);
    c_max    = R*a_today+z_today+b;
    %GSS and fminbnd return the *minimizer* of rhs_bellman
    obj      = @(c) -( u(c)+beta*myinterp1_equi(asset_grid,EVz,R*a_today+z_today-c) );
    optim    = GSS(obj,tiny,c_max-tiny,tol_gss);
    c_pol(ix) = optim;
    Tv(ix)    = -obj(optim);

end % end ix loop

c_pol = reshape(c_pol,[na,nz]);
Tv    = reshape(Tv,[na,nz]);

end %END FUNCTION <bellman_operator>

