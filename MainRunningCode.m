clc
close all
clear all
addpath(genpath('altmany-export_fig-83ee7fd\'))
addpath(genpath('denoisers\'))
%% Input parameters 
r            = im2double((imread('cameraman.tif')));
sigma_w      = 0.1;
A            = 'fft';
seedNum      = 1;
denoiserType = 'TV';
realOnly     = true;
figure, imshow(abs(r),[]), colorbar
title('Object-Reflectance')
set(gcf, 'Position', get(0, 'Screensize'));
export_fig('Figures\OpticalReflectance.png');
%% Generate w,g, and y form input-parameters
[M,N]    = size(r);                                                         % Size of the input-image 
g        = (sqrt(r)/2).*randn(M,N)+1j*(sqrt(r)/2).*randn(M,N);              % g ~ CN(0,D(r))
if(strcmp(A,'fft'))
    y = fftshift(fft2(g));
end
w        = (sqrt(sigma_w)/2)*randn([M,N])+1j*(sqrt(sigma_w)/2)*randn([M,N]);% w ~ CN(0,\sigma_w^2)   
y        = y+w;                                                             % noisy-measurements 
figure, imshow(log10(abs(y)),[]), colorbar
title('Noisy-Fourier domain (log)')
set(gcf, 'Position', get(0, 'Screensize'));
export_fig('Figures\NoisyFourierDomain.png');
figure, imshow(log10(abs(g)),[]), colorbar
title('Complex-optical field')
set(gcf, 'Position', get(0, 'Screensize'));
export_fig('Figures\ComplexOpticalField.png');
%% Plug and play ADMM algorithm     
maxIters     = 25;
v0           = abs(ifft2(y)).^2; 
u0           = zeros(size(v0));
sigmaLambda  = 0.5*sqrt(var(v0(:))); 
sigman       = 2; 
G            = denoiser(denoiserType,realOnly,sigman);
costFunction= zeros(maxIters,1);

vPrev        = v0;
uPrev        = u0; 
rPrev        = v0-u0;

for iters = 1:maxIters
    [c,mu]     = computeCovarianceAndMean(y,sigma_w,rPrev);
    rtildenext = vPrev-uPrev;
    figure(2),
    subplot(2,2,1), imshow(abs(rtildenext),[]), colorbar
    title('Input: Inversion-Op')
    costFunctionPrior = computeCostFunction(c,mu,rPrev,sigmaLambda,rtildenext);
    rNext             = inversionOperator(rtildenext,sigmaLambda,c,mu); 
    costFunctionPost(iters)  = computeCostFunction(c,mu,rNext,sigmaLambda,rtildenext);
    subplot(2,2,2), imshow(abs(rtildenext),[]), colorbar
    title('Output: Inversion-Op')
    vtildenext = rNext+uPrev;
    subplot(2,2,3), imshow(abs(vtildenext),[]), colorbar
    title('Input: Denoiser-Op')
    vNext      = G*vtildenext;
    subplot(2,2,4), imshow(abs(vNext),[]), colorbar
    title('Output: Denoiser-Op')
    
    uNext                       = uPrev+rNext-vNext;
    vPrev                       = vNext;
    uPrev                       = uNext;
    rPrev                       = rNext;
end
%%
figure, plot(costFunctionPost)
