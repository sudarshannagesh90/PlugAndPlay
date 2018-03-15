classdef denoiser < matlab.mixin.Copyable
% CLASS denoiser -     Denoise reconstructed objects. 
%              
% Constructor:
%   obj =
%   denoiser(denoiserType,sigmaAmplitude,sigmaPhase)
%
% Inputs:
%    denoiserType    : Enter the type of denoiser

% Output:
%    obj                     : denoiser type object.
%
% Overloaded Methods:
%    mtimes(obj,x)           : Called when used as (obj*x). Used to denoise the image with denoiser of the type denoiserType. 

% Author   : Sudarshan Nagesh             
% Institute: NorthWestern University (NU)

properties

    % Properties provided as input.
    denoiserType
    realOnly 
    sigman
    % Properties calculated by the class.

end

methods
    % Function constructFactorGraphMagicSquare - Constructor.
    function obj = denoiser(denoiserType,realOnly,sigman)
        obj.denoiserType        = denoiserType;                      
        obj.realOnly            = realOnly;
        obj.sigman              = sigman;
    end
    
    % Overloaded function for * (mtimes()).
    function [Result1] = mtimes(obj,x)
            AmplitudeEstimate = abs(x);
            if obj.realOnly == false
            PhaseEstimate     = angle(x);
            end
            
            switch obj.denoiserType
             case 'TV'
                options.lambda = obj.sigman;
                options.niter  = 25;
                options.display= 0;
                [RegularizedAmplitudeEstimate,~,~,~] = perform_tv_denoising(AmplitudeEstimate,options);
                if obj.realOnly == false
                    [RegularizedPhaseEstimate,~,~,~]     = perform_tv_denoising(PhaseEstimate,options);
                end
             case 'BM3D'
                [~, RegularizedAmplitudeEstimate] = BM3D(1, AmplitudeEstimate,25);
                if obj.realOnly == false
                    [~, RegularizedPhaseEstimate]     = BM3D(1, PhaseEstimate,25);
                end
             case 'DnCNN'
                noiseSigmaAmplitude              = 22; % noiseSigma must be an integer between 0 and 255. 
                [RegularizedAmplitudeEstimate]    = dnCNN(AmplitudeEstimate,noiseSigmaAmplitude);
                if obj.realOnly == false
                    noiseSigmaPhase                  = 25; % noiseSigma must be an integer between 0 and 255. 
                    [RegularizedPhaseEstimate]        = dnCNN(PhaseEstimate,noiseSigmaPhase);
                end
            end
         if obj.realOnly == false   
            Result1 = [RegularizedAmplitudeEstimate; RegularizedPhaseEstimate];
         else
            Result1 = [RegularizedAmplitudeEstimate;];
         end
    end   
end

end


