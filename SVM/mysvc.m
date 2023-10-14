function [nsv, alpha, b0] = mysvc(X,Y,ker,C,p1,p2)
%SVC Support Vector Classification
%
%  Usage: [nsv alpha bias] = svc(X,Y,ker,C)
%
%  Parameters: X      - Training inputs
%              Y      - Training targets
%              ker    - kernel function
%              C      - upper bound (non-separable case)
%              nsv    - number of support vectors
%              alpha  - Lagrange Multipliers
%              b0     - bias term
%
%  Author: Issam ElNaqa

% fprintf('Support Vector Classification\n')
% fprintf('_____________________________\n')
n = size(X,1);

% tolerance for Support Vector Detection
epsilon = svtol(C);

% Construct the Kernel matrix
%  H=y'*k(x,x)*y
YY=Y(:)*Y(:)';
H = YY.*mysvkernel(ker,X,X,p1,p2);

c = -ones(n,1);

% Add small amount of zero order regularisation to
% avoid problems when Hessian is badly conditioned.
H = H+1e-10*eye(size(H));

% Set up the parameters for the Optimisation problem
A=[];
b=[];
Aeq=[];
beq=[];

vlb = zeros(n,1);      % Set the bounds: alphas >= 0
vub = C*ones(n,1);     %                 alphas <= C
x0 = zeros(n,1);       % The starting point is [0 0 0   0]
neqcstr = nobias(ker); % Set the number of equality constraints (1 or 0)
if neqcstr
    Aeq = Y;, beq = 0;     % Set the constraint Ax = b
else
    A = [];, b = [];
end

% Solve the Optimisation Problem

% fprintf('Optimising ...\n');
st = cputime;

[alpha,FVAL,EXITFLAG,OUTPUT,lambda]=quadprog(H,c,A,b,Aeq,beq,vlb,vub,x0);
% fprintf('Execution time : %4.1f seconds\n',cputime - st);
% EXITFLAG
% OUTPUT

w2 = alpha'*H*alpha;
% fprintf('|w0|^2    : %f\n',w2);
% fprintf('Margin    : %f\n',2/sqrt(w2));
% fprintf('Sum alpha : %f\n',sum(alpha));


% Compute the number of Support Vectors
svi = find( alpha > epsilon);
nsv = length(svi);
% fprintf('Support Vectors : %d (%3.1f%%)\n',nsv,100*nsv/n);

% Implicit bias, b0
b0 = 0;
% check those changes...
% Explicit bias, b0
if nobias(ker) ~= 0
    % find b0 from average of support vectors on margin
    % SVs on margin have alphas: 0 < alpha < C
    svii = find( alpha > epsilon & alpha < (C - epsilon));
    if length(svii) > 0
        b0 =  (1/length(svii))*sum(Y(svii)' - H(svii,svi)*alpha(svi).*Y(svii)');
    else
        fprintf('No support vectors on margin - cannot compute bias.\n');
    end
end

end


