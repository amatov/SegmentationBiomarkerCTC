function [beta,bias, svi,nsv]=mysvrfnm(X1,Xo,Y1,Yo,ker,e,C1,C0)
%SVR Support Vector Regression
%  Parameters: X      - Training inputs
%              Y      - Training targets
%              ker    - kernel function
%              C      - upper bound (non-separable case)
%              loss   - loss function
%              e      - insensitivity
%              nsv    - number of support vectors
%              beta   - Difference of Lagrange Multipliers
%              bias   - bias term
% X,Y columns
% modified with a dual cost function...
%e=0.0;
loss='eInsensitive';
%C=Inf;
%ker='linear';
% tolerance for Support Vector Detection
epsilon1 = svtol(C1);
epsilono = svtol(C0);
epsilon=max(epsilon1, epsilono);
% construct the training data
X=[X1;Xo];Y=[Y1;Yo];n1=length(Y1);
no=length(Yo);
n=n1+no;
% Construct the Kernel matrix
H=mysvkernel(ker,X,X);
% Set up the parameters for the Optimisation proble
A=[]; 
b=[];
Aeq=[];  
beq=[];
switch lower(loss)
case 'einsensitive'
   Hb = [H -H; -H H];
   c = [(e*ones(n,1) - Y); (e*ones(n,1) + Y)];  
   vlb = zeros(2*n,1);    % Set the bounds: alphas >= 0
   vub1 = C1*ones(n1,1);  %       alphas <= C1 for class 1
   vubo = C0*ones(no,1);  %       alphas <= Co for class o   vub=[vub1;vubo];% alphas <= C
   x0 = zeros(2*n,1);     % The starting point is [0 0 0   0]
   neqcstr = nobias(ker); % Set the number of equality constraints (1 or 0)  
   if neqcstr
      Aeq = [ones(1,n) -ones(1,n)];, beq = 0;     % Set the constraint Ax = b
   else
      A = [];, b = [];
   end
case 'quadratic',
   Hb = H + eye(n)/(2*C);
   c = -Y;
   vlb = -1e30*ones(n,1);   
   vub = 1e30*ones(n,1);    
   x0 = zeros(n,1);              % The starting point is [0 0 0   0]
   neqcstr = nobias(ker);        % Set the number of equality constraints (1 or 0)  
   if neqcstr
      Aeq = ones(1,n);, beq = 0;      % Set the constraint Ax = b
   else
      A = [];, b = [];
   end
otherwise, disp('Error: Unknown Loss Function\n')
end

% Add small amount of zero order regularisation to 
% avoid problems when Hessian is badly conditioned. 
% Rank is always less than or equal to n.
% Note that adding to much reg will peturb solution

Hb = Hb+1e-10*eye(size(Hb));

% Solve the Optimisation Problem

fprintf('Optimising ...\n');
st = cputime;

%[alpha lambda how] = qp(Hb, c, A, b, vlb, vub, x0, neqcstr);
%options=optimset('LargeScale','off','Display','off');
% 'ShowStatusWindow','iterplus');
[alpha,FVAL,EXITFLAG,OUTPUT,lambda]=quadprog(Hb,c,A,b,Aeq,beq,vlb,vub,x0);
fprintf('Execution time : %4.1f seconds\n',cputime - st);
EXITFLAG
OUTPUT


switch lower(loss)
case 'einsensitive',
   beta =  alpha(1:n) - alpha(n+1:2*n);
case 'quadratic',
   beta = alpha;
end
fprintf('|w0|^2    : %f\n',beta'*H*beta);  
fprintf('Sum beta : %f\n',sum(beta));

% Compute the number of Support Vectors
svi = find( abs(beta) > epsilon );
nsv = length( svi );
fprintf('Support Vectors : %d (%3.1f%%)\n',nsv,100*nsv/n);

% Implicit bias, b0
bias = 0;

% Explicit bias, b0 
if nobias(ker) ~= 0
   switch lower(loss)
   case 'einsensitive',
      % find bias from average of support vectors with interpolation error e
      % SVs with interpolation error e have alphas: 0 < alpha < C
      svii = find( abs(beta) > epsilon & abs(beta) < (C - epsilon));
      if length(svii) > 0
         bias = (1/length(svii))*sum(Y(svii) - e*sign(beta(svii)) - H(svii,svi)*beta(svi));
      else 
         fprintf('No support vectors with interpolation error e - cannot compute bias.\n');
         bias = (max(Y)+min(Y))/2;
      end
   case 'quadratic',
      bias = mean(Y - H*beta);
   end 
end

% plotting the output
return