function test_anderson(m,deltas)
%test_anderson    Test Anderson acceleration for NCM problem.
%   test_anderson(M, DELTAS) computes the nearest correlation matrix with
%   smallest eigenvalue at least DELTAS(i) for two
%   invalid correlation test matrices by the unaccelerated alternating
%   projections method, nearcorr_new, and by alternating
%   projections with Anderson acceleration, nearcorr_aa, with 
%   history length M.
%   Vector DELTAS specifies different values of the lower bound
%   delta on the smallest eigenvalue of the solution.
%   Default: M = 2, DELTAS = [0 1e-8 0.1].

if nargin < 1, m = 2; end
if nargin < 2, deltas = [0 1e-8 0.1]; end

f = get(0,'Format'); % Save current format.
format shortg, format compact

fprintf('*** Test matrix 1:\n')
% Turkay, Epperlein, and Christofides (2003).
A = [1 -0.55 -0.15 -0.10;
    -0.55 1 0.90 0.90;
    -0.15 0.90 1 0.90;
    -0.10 0.90 0.90 1]

pattern = [];

for delta = deltas
    fprintf('--> m = %d, delta = %d\n', m, delta)
    [X,it] = nearcorr_new(A,pattern,delta);
    [XAA,itAA] = nearcorr_aa(A,pattern,m,100,'u',delta); 
    fprintf(['Iterations for NCM: standard = %d, accelerated = %d.' ...
             '  Reduction factor = %.2f\n'], it, itAA, it/itAA)

end

fprintf('\n*** Test matrix 2 and the pattern of elements to keep fixed:\n')
% Finger (1997).
A = [1 0.18 -0.13 -0.26 0.19 -0.25 -0.12;
    0.18 1 0.22 -0.14 0.31 0.16 0.09;
   -0.13 0.22 1 0.06 -0.08 0.04 0.04;
   -0.26 -0.14 0.06 1 0.85 0.85 0.85;
    0.19 0.31 -0.08 0.85 1 0.85 0.85;
   -0.25 0.16 0.04 0.85 0.85 1 0.85;
   -0.12 0.09 0.04 0.85 0.85 0.85 1]

% Define pattern to keep fixed
pattern = blkdiag(ones(3,3),eye(4))

for delta = deltas
    fprintf('--> m = %d, delta = %d\n', m, delta)
    [X,it] = nearcorr_new(A,pattern,delta);
    [XAA,itAA] = nearcorr_aa(A,pattern,m,100,'u',delta); 
    fprintf(['Iterations for NCM: standard = %d, accelerated = %d.' ...
             '  Reduction factor = %.2f\n'], it, itAA, it/itAA)

end

format(f) % Restore original format.
