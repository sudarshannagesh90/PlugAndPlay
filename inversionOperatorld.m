function [rnext,costRoot] = inversionOperatorld(rtilde,sigmaLambda,c,mu)
[M,N] = size(rtilde);
rtilde= rtilde(:);
rnext = zeros(size(rtilde));

for index = 1:size(c,1)*size(c,2)
alpha1             = 1/sigmaLambda^2;
alpha2             = -rtilde(index)/sigmaLambda^2;
alpha3             = 1;
alpha4             = -(c(index)+mu(index)^2);
rootsEquation      = roots([alpha1 alpha2 alpha3 alpha4]);

temp=rootsEquation;
temp_diffs=abs(temp-real(temp));
temp_i=temp(temp_diffs==min(temp_diffs));
root=real(temp_i(1));
costRoot(index) = computeCostFunction(c(index),mu(index),root,sigmaLambda,rtilde(index));

for ind = 1:length(rootsEquation)
rootsCost(ind)     = computeCostFunction(c(index),mu(index),rootsEquation(ind),sigmaLambda,rtilde(index));
end
[~,minPos]         = min(rootsCost);
rootsEquation      = abs(rootsEquation(minPos));
rnext(index)       = rootsEquation;
end
rnext = reshape(rnext,M,N);
