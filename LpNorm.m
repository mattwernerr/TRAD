function fNormed = LpNorm(x, f, p)
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% Calculate the p-norm of the function f(x) by the standard definition
% 
%                /   _            \  1/p
%               |   /              |
%     ||f||  =  |   |        p     |
%          p    |   |  |f(x)|  dx  |
%                \ _/             /
% 
% 
% 
%    Inputs:
% 
%                 x - Independent variable where f is evaluated.
%                     Size: n-by-1 (vector)
% 
%                 f - Dependent variable evaluated at x.
%                     Size: n-by-1 (vector)
% 
%                 p - Positive number, including INF, indicating which
%                     norm to be taken.
%                     Size: 1-by-1 (scalar)
% 
%    Outputs:
% 
%           fnormed - Calculated number of f's p-norm.
%                     Size: 1-by-1 (scalar)
% 

% Check that p is positive
assert(p > 0, "p must be positive.")

% Handle the p = infinity case
if (isinf(p))
    fNormed = max(f);
    return
end

% Calculate for all other p
fNormed = trapz(x, abs(f).^p).^(1/p);