function  logLikelihood = computeLoglikelihoodFunction(c,mu,r)

epsilon = 1e-5;
logLikelihood = sum(log10((r(:)+epsilon)))+sum(sum((1./r(:)).*(c(:)+abs(mu(:)).^2)));