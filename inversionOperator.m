function [rnext] = inversionOperator(rtilde,sigmaLambda,c,mu)
rnext = zeros(size(rtilde));
for index = 1:size(c,1)*size(c,2)
alpha1             = 1/sigmaLambda^2;
alpha2             = -rtilde(index)/sigmaLambda^2;
alpha3             = 1;
alpha4             = -(c(index)+mu(index));
rootsEquation      = roots([alpha1 alpha2 alpha3 alpha4]);
rootsEquation      = min(abs(rootsEquation));
rnext(index)       = abs(rootsEquation);
end
