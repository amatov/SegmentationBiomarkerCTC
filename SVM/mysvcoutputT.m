function predictedY = mysvcoutputT(trnX,trnY,tstX,ker,p1,p2,alpha,bias)
%SVCOUTPUT Calculate SVC Output
%
%  Usage: predictedY = svcoutput(trnX,trnY,tstX,ker,alpha,bias,actfunc)
%
%  Parameters: trnX   - Training inputs
%              trnY   - Training targets
%              tstX   - Test inputs
%              ker    - kernel function
%              beta   - Lagrange Multipliers
%              bias   - bias              
%              actfunc- activation function (0(default) hard | 1 soft) 
%              
%              No thresholding
%  Author: Issam ElNaqa

if (nargin < 6) % check correct number of arguments
   help svcoutputT
else
   %ker='rbf';
   trnYY=ones(size(tstX,1),1)*trnY(:)';
   H=trnYY .* mysvkernel(ker,tstX,trnX,p1,p2);
   %H = svkernelclassout(ker,tstX,trnX,trnY,p);
   predictedY=H*alpha+bias;
end