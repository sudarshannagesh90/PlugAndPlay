function [rnext,costFunction] = inversionOperator(rtilde,sigmaLambda,c,mu)
[M,N] = size(rtilde);
rtilde= rtilde(:);
rnext = zeros(size(rtilde));
costFunction = zeros(size(rtilde));
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
costRoot = computeCostFunction(c(index),mu(index),root,sigmaLambda,rtilde(index));
if costRoot>100
x = -10:0.01:10;
for ind = 1:length(x)
    polynomialRoots(ind) = alpha1*x(ind)^3+alpha2*x(ind)^2+alpha3*x(ind)+alpha4;
end
figure(100), plot(x,polynomialRoots)
hold on
plot(x,zeros(size(x)),'r--')
title(['Cost: ',num2str(costRoot)])
pause
end
for ind = 1:length(rootsEquation)
rootsCost(ind)     = computeCostFunction(c(index),mu(index),rootsEquation(ind),sigmaLambda,rtilde(index));
end
[costFunction(index),minPos]   = min(rootsCost);
rootsEquation      = abs(rootsEquation(minPos));
rnext(index)       = rootsEquation;
end
rnext = reshape(rnext,M,N);
costFunction = reshape(costFunction,M,N);
