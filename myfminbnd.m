function [xf,fval,exitflag] = myfminbnd(funfcn,ax,bx,tol,maxfun,maxiter,varargin)
%FMINBND Single-variable bounded nonlinear function minimization.
%   X = FMINBND(FUN,x1,x2) attempts to find  a local minimizer X of the function
%   FUN in the interval x1 < X < x2.  FUN is a function handle.  FUN accepts
%   scalar input X and returns a scalar function value F evaluated at X.
%
%   X = FMINBND(FUN,x1,x2,OPTIONS) minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, created with
%   the OPTIMSET function. See OPTIMSET for details. FMINBND uses these
%   options: Display, TolX, MaxFunEval, MaxIter, FunValCheck, PlotFcns,
%   and OutputFcn.
%
%   X = FMINBND(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with the function FUN in PROBLEM.objective, the interval
%   in PROBLEM.x1 and PROBLEM.x2, the options structure in PROBLEM.options,
%   and solver name 'fminbnd' in PROBLEM.solver.
%
%   [X,FVAL] = FMINBND(...) also returns the value of the objective function,
%   FVAL, computed in FUN, at X.
%
%   [X,FVAL,EXITFLAG] = FMINBND(...) also returns an EXITFLAG that
%   describes the exit condition. Possible values of EXITFLAG and the
%   corresponding exit conditions are
%
%    1  FMINBND converged with a solution X based on OPTIONS.TolX.
%    0  Maximum number of function evaluations or iterations reached.
%   -1  Algorithm terminated by the output function.
%   -2  Bounds are inconsistent (that is, ax > bx).
%
%   [X,FVAL,EXITFLAG,OUTPUT] = FMINBND(...) also returns a structure
%   OUTPUT with the number of iterations taken in OUTPUT.iterations, the
%   number of function evaluations in OUTPUT.funcCount, the algorithm name
%   in OUTPUT.algorithm, and the exit message in OUTPUT.message.
%
%   Examples
%     FUN can be specified using @:
%        X = fminbnd(@cos,3,4)
%      computes pi to a few decimal places and gives a message upon termination.
%        [X,FVAL,EXITFLAG] = fminbnd(@cos,3,4,optimset('TolX',1e-12,'Display','off'))
%      computes pi to about 12 decimal places, suppresses output, returns the
%      function value at x, and returns an EXITFLAG of 1.
%
%     FUN can be an anonymous function:
%        X = fminbnd(@(x) sin(x)+3,2,5)
%
%     FUN can be a parameterized function.  Use an anonymous function to
%     capture the problem-dependent parameters:
%        f = @(x,c) (x-c).^2;  % The parameterized function.
%        c = 1.5;              % The parameter.
%        X = fminbnd(@(x) f(x,c),0,1)
%
%   See also OPTIMSET, FMINSEARCH, FZERO, FUNCTION_HANDLE.

%   References:
%   "Algorithms for Minimization Without Derivatives",
%   R. P. Brent, Prentice-Hall, 1973, Dover, 2002.
%
%   "Computer Methods for Mathematical Computations",
%   Forsythe, Malcolm, and Moler, Prentice-Hall, 1976.

%   Original coding by Duane Hanselman, University of Maine.
%   Copyright 1984-2019 The MathWorks, Inc.

% Set default options
% options.TolX = 1e-8;
% options.MaxFunEvals  = 500;
% options.MaxIter = 500;
% options.printtype = 'notify';

%tol       = options.TolX;
%maxfun    = 500;
%maxiter   = 500;
%print     = 1;  % possible values 0,1,2,3 see switch below

funccount = 0;
iter      = 0;

% switch printtype
%     case {'notify','notify-detailed'}
%         print = 1;
%     case {'none','off'}
%         print = 0;
%     case {'iter','iter-detailed'}
%         print = 3;
%     case {'final','final-detailed'}
%         print = 2;
%     otherwise
%         print = 1;
% end

% checkbounds
if ax > bx
    exitflag = -2;
    xf=[]; fval = [];
    msg='Exiting: Lower Bound Exceeds Upper Bound';
    %output.iterations = 0;
    %output.funcCount = 0;
    %output.algorithm = 'golden section search, parabolic interpolation';
    %output.message = msg;
    warning(msg)
    % Have not initialized OutputFcn; do not need to call it before returning
    return
end

% Assume we'll converge
exitflag = 1;

%header = ' Func-count     x          f(x)         Procedure';
%procedure='       initial';

% Compute the start point
seps = sqrt(eps);
c = 0.5*(3.0 - sqrt(5.0));  % 0.381966011250105

a = ax; b = bx;
v = a + c*(b-a);
w = v; 
xf = v;
d = 0.0; e = 0.0;
x= xf; 
fx = funfcn(x,varargin{:});
funccount = funccount + 1;

% Check that the objective value is a scalar
if numel(fx) ~= 1
    error('The objective function is not scalar');
end

% Display the start point if required
%if print > 2
%    disp(' ')
%    disp(header)
%    fprintf('%5.0f   %12.6g %12.6g %s\n',funccount,xf,fx,procedure)
%end

fv = fx; fw = fx;
xm = 0.5*(a+b);
tol1 = seps*abs(xf) + tol/3.0;
tol2 = 2.0*tol1;

% Main loop
while ( abs(xf-xm) > (tol2 - 0.5*(b-a)) )
    gs = 1;
    % Is a parabolic fit possible
    if abs(e) > tol1
        % Yes, so fit parabola
        gs = 0;
        r = (xf-w)*(fx-fv);
        q = (xf-v)*(fx-fw);
        p = (xf-v)*q-(xf-w)*r;
        q = 2.0*(q-r);
        if q > 0.0,  p = -p; end
        q = abs(q);
        r = e;  e = d;
        
        % Is the parabola acceptable
        if ( (abs(p)<abs(0.5*q*r)) && (p>q*(a-xf)) && (p<q*(b-xf)) )
            
            % Yes, parabolic interpolation step
            d = p/q;
            x = xf+d;
            %procedure = '       parabolic';
            
            % f must not be evaluated too close to ax or bx
            if ((x-a) < tol2) || ((b-x) < tol2)
                si = sign(xm-xf) + ((xm-xf) == 0);
                d = tol1*si;
            end
        else
            % Not acceptable, must do a golden section step
            gs=1;
        end
    end
    if gs
        % A golden-section step is required
        if xf >= xm
            e = a-xf;
        else
            e = b-xf;
        end
        d = c*e;
        %procedure = '       golden';
    end
    
    % The function must not be evaluated too close to xf
    si = sign(d) + (d == 0);
    x = xf + si * max( abs(d), tol1 );
    fu = funfcn(x,varargin{:});
    funccount = funccount + 1;
    
    iter = iter + 1;
    %if print > 2
    %    fprintf('%5.0f   %12.6g %12.6g %s\n',funccount, x, fu, procedure);
    %end
    
    % Update a, b, v, w, x, xm, tol1, tol2
    if fu <= fx
        if x >= xf
            a = xf;
        else
            b = xf;
        end
        v = w; fv = fw;
        w = xf; fw = fx;
        xf = x; fx = fu;
    else % fu > fx
        if x < xf
            a = x;
        else
            b = x;
        end
        if ( (fu <= fw) || (w == xf) )
            v = w; fv = fw;
            w = x; fw = fu;
        elseif ( (fu <= fv) || (v == xf) || (v == w) )
            v = x; fv = fu;
        end
    end
    xm = 0.5*(a+b);
    tol1 = seps*abs(xf) + tol/3.0; tol2 = 2.0*tol1;
    
    if funccount >= maxfun || iter >= maxiter
        exitflag = 0;
        %output.iterations = iter;
        %output.funcCount = funccount;
        %output.algorithm = 'golden section search, parabolic interpolation';
        fval = fx;
        return
    end
end % while

fval = fx;
%output.iterations = iter;
%output.funcCount = funccount;
%output.algorithm = 'golden section search, parabolic interpolation';

end %end function

