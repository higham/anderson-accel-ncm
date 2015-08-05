function [X,iter] = nearcorr_aa(A,pattern,mMax,itmax,ls_solve,delta,...
                                tol,droptol,beta,AAstart)
%nearcorr_aa     Nearest correlation matrix with Anderson Acceleration.
%   [X,ITER] = nearcorr_aa(A,PATTERN,MMAX,ITMAX,LS_SOLVE,DELTA,...
%                          TOL,DROPTOL,BETA,AASTART)
%   finds the nearest correlation matrix X to a symmetric matrix A by
%   the alternating projections method with Anderson acceleration.
%   ITER is the number of iterations taken.
%   Prescribed elements can be kept fixed by specifying them in the
%   matrix PATTERN: the (i,j) element is 1 if the corresponding element 
%   of A is to remain fixed, else zero. 
%   By default PATTERN is empty, meaning no elements are fixed.
%   If PATTERN is non-empty then the unit diagonal must be explicitly
%   forced if required.
%   The solution can forced to be positive definite with smallest
%   eigenvalue at least DELTA, 0 < DELTA <= 1. Default: DELTA = 0.
%   MMAX = history length parameter (non-negative integer).
%          Default: 2.  MMAX = 0 means no acceleration.
%   ITMAX = maximum allowable number of iterations.  Default: 100.
%   LS_SOLVE = the least-squares solve method:
%        'u' -- QR factorization with updating (original code, default),
%        'n' -- normal equations,
%        'b' -- MATLAB's backslash.
%   TOL = convergence tolerance. Default: length(X)*eps.
%   DROPTOL = tolerance for dropping stored residual vectors to improve
%           conditioning: If DROPTOL > 0, drop residuals if the
%           condition number exceeds DROPTOL; if DROPTOL <= 0,
%           do not drop residuals. Default: 0.
%   BETA = damping factor: If BETA > 0 (and BETA ~= 1), then the step is
%       damped by BETA; otherwise, the step is not damped (default).
%       NOTE: BETA can be a function handle; form beta(iter), where iter
%       is the iteration number and 0 < beta(iter) <= 1.
%   AASTART = acceleration delay factor: If AASTART > 0, start acceleration
%         when iter = AAstart. Default: 1.

%   Based on the code AndAcc.m by Homer Walker dated 10/14/2011,
%   which is listed in 
%   Homer F. Walker, Anderson acceleration: Algorithms and implementations.
%   Technical Report MS-6-15-50, Mathematical Sciences Department, 
%   Worcester Polytechnic Institute, Worcester, MA 01609, USA, June 2011.
%   The code is specialized for the alternating projections method for the
%   nearest correlation matrix, adds two extra options for solving
%   the least squares problem, and contains many other changes.

%   Nick Higham and Natasa Strabic, 2015.

if ~isequal(A,A'), error('The matrix A must by symmetric.'), end
n = length(A);

% Set the method parameters.
if nargin < 2, pattern = []; end
if nargin < 3 || isempty(mMax), mMax = 2; end
if nargin < 4 || isempty(itmax), itmax = 100; end
if nargin < 5 || isempty(ls_solve), ls_solve = 'u'; end
if nargin < 6 || isempty(delta), delta = 0; end
if nargin < 7 || isempty(tol), tol = n*eps; end
if nargin < 8 || isempty(droptol), droptol = 0; end
if nargin < 9 || isempty(beta), beta = 1; end
if nargin < 10, AAstart = 1; end

% Initialize the storage arrays.
DG = []; % Storage of g-value differences.
DF = []; % Storage of f-value differences, needed for 'n' and 'b' options.

% Initialize the number of stored residuals.
mAA = 0;

% Initialization for the nearest corr. matrix problem.
Yin = A;
Sin = zeros(n);

% Top of the iteration loop.
for iter = 1:itmax
    % Apply g
    [Xout,Yout,Sout] = ap_step(A,Yin,Sin,pattern,delta);
    
    % We have to vectorize the output
    x = [Yin(:); Sin(:)];
    gval = [Yout(:); Sout(:)];

    fval = gval - x;
    
    % Convergence test
    rel_diffXY = norm(Yout - Xout,'fro')/norm(Yout,'fro');
    if rel_diffXY < tol, break; end
    
    if mMax == 0 || iter < AAstart,
    % Without acceleration, update x <- g(x) to obtain the next
    % approximate solution.
        x = gval;
    else
    % Apply Anderson acceleration.
    
    % Update the df vector and the DG array.
    % If the least-squares solves are by normal equations or backslash, we
    % need DF, the matrix of function value differences built in the same
    % way as DG.
        if iter > AAstart,
            df = fval-f_old;
            if mAA < mMax,
                DG = [DG gval-g_old];
                if ls_solve ~= 'u'
                    DF = [DF df];
                end
            else
                DG(:,1)=[];
                DG(:,mAA) = gval-g_old;
                if ls_solve ~= 'u'
                    DF(:,1) = [];
                    DF(:,mAA) = df;
                end
            end
            mAA = mAA + 1;
        end
        f_old = fval;
        g_old = gval;
        
        if mAA == 0
            % If mAA == 0, update x <- g(x) to obtain the next approximate 
            % solution.  This is the x1 = g(x0) line in pseudocode.
            x = gval;
        else
        % If mAA > 0, solve the least-squares problem and update the
        % solution.        
            switch ls_solve
              case 'u' % Use the QR factors.
                if mAA == 1
                % If mAA == 1, form the initial QR decomposition.
                    R(1,1) = norm(df);
                    Q = R(1,1)\df;
                else
                % If mAA > 1, update the QR decomposition.
                    if mAA > mMax
                    % If the column dimension of Q is mMax, delete
                    % the first column and update the decomposition.
                        [Q,R] = qrdelete(Q,R,1);
                        mAA = mAA - 1;
                    % The following treats the qrdelete quirk described below.
                        if size(R,1) ~= size(R,2),
                            Q = Q(:,1:mAA-1); R = R(1:mAA-1,:);
                        end
                    % Explanation: If Q is not square, then
                    % qrdelete(Q,R,1) reduces the column
                    % dimension of Q by 1 and the column and row
                    % dimensions of R by 1. But if Q *is* square, then
                    % the column dimension of Q is not reduced and only
                    % the column dimension of R is reduced by one. This
                    % is to allow for MATLAB’s default "thick" QR
                    % decomposition, which always produces a square Q.
                    end
                    % Now update the QR decomposition to incorporate
                    % the new column.
                    for j = 1:mAA - 1
                        R(j,mAA) = Q(:,j)'*df;
                        df = df - R(j,mAA)*Q(:,j);
                    end
                    R(mAA,mAA) = norm(df);
                    Q = [Q,R(mAA,mAA)\df];
                end

                if droptol > 0
                % Drop residuals to improve conditioning if necessary.
                    condDF = cond(R);
                    while condDF > droptol && mAA > 1
                        fprintf('cond(D) = %e, reducing mAA to %d \n', ...
                                condDF, mAA-1);
                        [Q,R] = qrdelete(Q,R,1);
                        DG = DG(:,2:mAA);
                        mAA = mAA - 1;
                        % Treat the qrdelete quirk described above.
                        if size(R,1) ~= size(R,2),
                            Q = Q(:,1:mAA); R = R(1:mAA,:);
                        end
                        condDF = cond(R);
                    end
                end % if droptolerance

                % Solve the least-squares problem.
                gamma = R\(Q'*fval);

              case 'n' % use normal equations
                  if mAA > mMax, mAA = mAA - 1; end
                  c = DF'*fval;
                  [~,R] = qr(DF,0); % economy size QR
                  y = R'\c;
                  gamma = R\y;

              case 'b' % use backslash
                  if mAA > mMax, mAA = mAA - 1; end
                  gamma = DF\fval;
            end % switch                   
                                      
            % Update the approximate solution.
            x = gval - DG*gamma;            
  
        end % end of least-squares solve, we have the next iterate x
        
    end % end of the acceleration part
    
    % Update the parameters: undo reshape of x to get new Yin, Sin.
    Yin(:) = x(1:n^2);
    Sin(:) = x(n*n+1:2*n^2);

end
% Bottom of the iteration loop.

if rel_diffXY > tol && iter == itmax,
    error(['Stopped after ' num2str(itmax) ' its. Try increasing ITMAX.'])
end

X = Yin;

end

% --------------------------------------------
% Subfunctions
function [Xout,Yout,Sout] = ap_step(A,Yin,Sin,pattern,delta)
% Fixed point iteration for the alternating projections step
% for the nearest correlation matrix.

R = Yin - Sin;
Xout = proj_spd(R,delta);
Sout = Xout - R;
Yout = proj_pattern(A,Xout,pattern);
end

function X = proj_pattern(A,X,pattern)
% Return the nearest matrix to X with fixed elements from A specified in
% the positions PATTERN
 
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
% Return the nearest positive semidefinite matrix to A with the smallest
% eigenvalue at least delta
 
[V,D] = eig(A);
X = V*diag(max(diag(D),delta))*V';
X = (X+X')/2; % Ensure symmetry.
end