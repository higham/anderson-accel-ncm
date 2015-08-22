function [Y,iter] = nearcorr_new(A,pattern,delta,tol,itmax)
% nearcorr_new  Nearest correlation matrix by alternating projections.
%   [Y,ITER] = nearcorr_new(A,PATTERN,DELTA,TOL,ITMAX) finds the
%   nearest correlation matrix Y to a symmetric matrix A by the 
%   alternating projections method.  ITER is the number of iterations taken.
%   Prescribed elements can be kept fixed by specifying them in the
%   matrix PATTERN: the (i,j) element is 1 if the corresponding element 
%   of A is to remain fixed, else zero.
%   By default PATTERN is empty, meaning no elements are fixed.
%   If PATTERN is non-empty then the unit diagonal must be explicitly
%   forced if required.
%   The solution can be forced to be positive definite with smallest
%   eigenvalue at least DELTA, 0 < DELTA <= 1. By default, DELTA = 0.
%   TOL is the convergence tolerance. Default: length(A)*eps.
%   ITMAX is the maximum number of iterations allowed. Default: 100.

%   This is a modified version of nearcorr.m by Nick Higham, from
%   https://nickhigham.wordpress.com/2013/02/13/the-nearest-correlation-matrix/
%   This version removes the support for weights and "highly
%   non-positive definite A" but adds the PATTERN and DELTA options.

%   Nick Higham and Natasa Strabic, 2015.

if ~isequal(A,A'), error('The input matrix must by symmetric.'), end
if nargin < 2, pattern = []; end
if nargin < 3 || isempty(delta), delta = 0; end
if nargin < 4 || isempty(tol), tol = length(A)*eps; end
if nargin < 5, itmax = 100; end

Y = A;
iter = 0;
rel_diffXY = inf;
dS = zeros(size(A));

while rel_diffXY > tol

    [X,Y,dS] = ap_step(A,Y,dS,pattern,delta);
    rel_diffXY = norm(Y-X,'fro')/norm(Y,'fro');
   
    iter = iter + 1;
    if iter > itmax
       error(['Stopped after ' num2str(itmax) ' its. Try increasing ITMAX.'])
    end
    
end
end

% --------------------------------------------
% Subfunctions

function [X,Y,dS] = ap_step(A,Y,dS,pattern,delta)
% Fixed point iteration for the alternating projections step 
% for the nearest correlation matrix.

R = Y - dS;
X = proj_spd(R,delta);
dS = X - R;
Y = proj_pattern(A,X,pattern);
end

function X = proj_pattern(A,X,pattern)
% Nearest matrix to X with fixed elements from A specified in
% the positions PATTERN.
 
n = length(A);
if isempty(pattern)
    % Standard NCM, with no fixed elements.  Set unit diagonal.
    X(1:n+1:n^2) = 1;
else
    F = find(pattern); % Locations of fixed elements in A.
    X(F) = A(F);
end
end

function X = proj_spd(A,delta)
% Return the nearest positive semidefinite matrix to A with thesmallest
% eigenvalue at least delta.
 
[V,D] = eig(A);
X = V*diag(max(diag(D),delta))*V';
X = (X+X')/2; % Ensure symmetry.
end
