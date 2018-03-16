function  costFunction = computeCostFunction(c,mu,r,sigmaLambda,rtilde)

epsilon = 1e-5;
costFunction = sum(log10(r(:)+epsilon))+sum(sum((1./r(:)).*(c(:)+abs(mu(:)).^2)))+(1/(2*sigmaLambda^2))*norm(r-rtilde)^2;
