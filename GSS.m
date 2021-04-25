function min = GSS(fx,a,b,tol)
% Golden Section Search
% fx: f(x) with a single variable defined by anonymous function
% a,b: [a,b] initial interval where f(x) is quasiconvex

% min: global minimizer of f(x) on [a,b]
r = (3-sqrt(5))/2; 
c = a + r*(b-a);
d = a + (1-r)*(b-a);

fc = fx(c);
fd = fx(d);

distance = 100;
while distance > tol^2
    if fc >= fd  % [a c d b] => [a c) is eliminated 
        x = c + (1-r)*(b-c);
        a = c;
        c = d; fc = fd;
        d = x; fd = fx(d);
    else% [a c d b] => (d b] is eliminated
        x = a + r*(d-a);
        b = d;
        d = c; fd = fc;
        c = x; fc = fx(c);
    end
    
    distance = d-c;
end
min = (c+d)/2;

