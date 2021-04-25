function g_star = compute_fixed_point(g,tol,max_iter,damp,method,verbose,par)
% This function computes the fixed point of g'=F(g)
% where G(.) is a functional operator. It can be either the
% Bellman operator or the Coleman operator
% Author: Alessandro Di Nola

%{
%INPUTS
g:        initial condition (a suitable function)
tol:      tolerance criterion
max_iter: maximum number of iterations allowed
damp:     dampening (weight to old iterate)
method:   string with method, bellman, coleman, coleman_python, fpi
verbose:  flag 0-1 controls verbosity
par:      matlab struct with parameters to be passed to the operator
%}

iter = 0;
dist = tol+1;

tic
disp(['Compute fixed point of ' method]);
fprintf(' \n')
fprintf('Starting Iterations... \n');
while dist>tol && iter<=max_iter
    iter = iter+1;
    
    switch method
        case 'bellman'
            %Apply the Bellman's operator
            T_g = bellman_operator(g,par);
            
        case 'coleman'
            %Apply the Coleman operator
            T_g = coleman_operator(g,par);
            
        case 'coleman_python'
            %Apply Coleman operator as rootfinding
            T_g = coleman_operator_python(g,par);
        case 'fpi'
            %Apply method suggested by Maffezzoli
            T_g = fixed_point_iteration(g,par);
            
        otherwise
            error('Method selected does not exist!')
    end
    
    %Compute sup-norm b/w two successive iterations
    dist = max(abs(T_g(:)-g(:)));
    %update current guess
    g = (1-damp)*T_g+damp*g; 
    %Display some output
    if verbose==1; fprintf('iter: %d, dist: %f \n',iter,dist); end
    
end %end while

if strcmp(method,'bellman')
    fprintf('Value Function Iteration completed after %d iterations \n',iter);
elseif strcmp(method,'coleman') || strcmp(method,'coleman_python')
    fprintf('Time Iteration completed after %d iterations \n',iter);
elseif strcmp(method,'fpi')
    fprintf('Fixed point iteration completed after %d iterations \n',iter);
end

if iter>=max_iter
    warning('Maximum number of iterations reached!')
end

toc

g_star = T_g; %can be either value or policy function depending on the method

end %END function <compute_fixed_point>
