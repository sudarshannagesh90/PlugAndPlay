classdef ForwardModelGenerator < matlab.mixin.Copyable
% CLASS FORWARDMODELGENERATOR -  Generates the forward-model coherent
%                                images under the assumption of the given 
%                                forward-model. 
%                                  .
%              
% Constructor:
%   obj = ForwardModelGenerator(objectSizePixels,measurementType,forwardModelType,seedNum,noiseType,sigma_w)
%
% Inputs:
%    objectSizePixels          : Size of the object in pixels given as [rows,cols]. 
%    measurementType           : Type of measurements: Can be one of linear or phase-less.
%    forwardModelType          : Type of forward-model used: Can be A=I(for identity) or A=F (for Fourier).  
%    seedNum                   : Sets the seed of the simulation (should be an integer). 
%    noiseType                 : Noise type of the measurement (can be one of 'Poisson' or 'Gaussian').
%    sigma_w                   : Std. deviation of the additive white-Gaussian noise. (optional)
%
% Output:
%    obj                       : ForwardModelGenerator type object with 
%                              overloaded functions.
%
% Overloaded Methods:
%    conj(obj)      : Returns a copy of the object.
%    transpose(obj) : Returns a copy of the operator.
%    ctranspose(obj): Returns a copy of the operator.
%    mtimes(obj,x)  : Called when used as (obj*x). Applies the specified 
%                     ForwardModelGenerator on an image x after reshaping it to
%                     objectSizePixels. The number of elements in x must be
%                     prod(objectSizePixels).
%

% Author   : Sudarshan Nagesh             
% Institute: NorthWestern University (NU) 
properties(Constant)
    listOfSupportedmeasurementType  = {'Linear','PhaseLess'};
    listOfSupportedforwardModelType = {'Identity','Fourier'};
    listOfSupportednoiseType        = {'Poisson','Gaussian'};
end


properties

    % Properties provided as input.
    objectSizePixels
    measurementType
    forwardModelType
    seedNum
    noiseType
    sigma_w
end

methods
    % Function detectorSamplingOperator - Constructor.
    function obj = ForwardModelGenerator(objectSizePixels,measurementType,forwardModelType,seedNum,noiseType,sigma_w)

        obj.objectSizePixels                       = objectSizePixels;
        obj.measurementType                        = measurementType;
        obj.forwardModelType                       = forwardModelType;
        obj.seedNum                                = seedNum;
        obj.noiseType                              = noiseType;
        obj.sigma_w                                = sigma_w;

    end
    
    % Overloaded function for conj().
    function res = conj(obj)
        res = obj;
    end
    
    % Overloaded function for .' (transpose()).
    function res = transpose(obj)
         res = obj;
    end
    
    % Overloaded function for ' (ctranspose()).
    function res = ctranspose(obj)
        res = obj;
    end
    
    % Overloaded function for * (mtimes()).
    function res = mtimes(obj,x)
        rng(obj.seedNum)
        x   = reshape(x,obj.objectSizePixels);
        x   = double(x);
        g        = sqrt(x/2).*randn(obj.objectSizePixels(1),obj.objectSizePixels(2))+1j*sqrt(x/2).*randn(obj.objectSizePixels(1),obj.objectSizePixels(2));   % g ~ CN(0,D(x))
        if(strcmp(obj.measurementType,'Linear'))
            if(strcmp(obj.forwardModelType,'Fourier'))
                y = fftshift(fft2(g));
            elseif (strcmp(obj.forwardModelType,'Identity'))
                y = g;
            end
        end
        if (strcmp(obj.noiseType,'Gaussian'))
            w        = sqrt(obj.sigma_w/2)*randn([obj.objectSizePixels(1),obj.objectSizePixels(2)])+1j*sqrt(obj.sigma_w/2)*randn([obj.objectSizePixels(1),obj.objectSizePixels(2)]);    % w ~ CN(0,\sigma_w^2)
            y        = y+w;  % noisy-measurements 
        elseif (strcmp(obj.noiseType,'Poisson'))
            y        = poissrnd(y);
        end
        res = [y(:);g(:)];
    end
end

methods    
    % Check validity of properties provided as input.
    
    function set.objectSizePixels(obj,objectSizePixels)
        validateattributes(objectSizePixels,...
                           {'numeric'},...
                           {'nonsparse','vector','numel',2,'integer','positive'},...
                           mfilename,'objectSizePixels',1);
        if ~isa(objectSizePixels,'double')
            objectSizePixels = double(objectSizePixels);
        end
        if ~isrow(objectSizePixels)
            objectSizePixels = objectSizePixels(:)';
        end
        obj.objectSizePixels = objectSizePixels;
    end
    
    function set.seedNum(obj,seedNum)
        validateattributes(seedNum,...
                           {'double','single'},...
                           {'nonsparse','scalar','real','nonnan','finite','positive',},...
                           mfilename,'seedNum');
        if ~isa(seedNum,'double')
            seedNum = double(seedNum);
        end
        obj.seedNum = seedNum;
    end
    function set.sigma_w(obj,sigma_w)
        validateattributes(sigma_w,...
                           {'double','single'},...
                           {'nonsparse','scalar','real','nonnan','finite','positive',},...
                           mfilename,'sigma_w');
        if ~isa(sigma_w,'double')
            sigma_w = double(sigma_w);
        end
        obj.sigma_w = sigma_w;
    end
    function set.measurementType(obj,measurementType)
        if ~isempty(measurementType)
            validateattributes(measurementType,...
                               {'char'},{'scalartext','nonempty'},...
                               mfilename,'transform',1);
            if ~ismember(measurementType,obj.listOfSupportedmeasurementType)
                error(strcat('Variable measurementType contains a method that is not supported.\n',...
                             'Supported measurement types are: %s.'),...
                             strjoin(cellfun(@(x) sprintf('''%s''',x),obj.listOfSupportedmeasurementType,'UniformOutput',false),', '));
            end
        else
            measurementType = 'Linear';
        end
        obj.measurementType = measurementType;
    end
    function set.forwardModelType(obj,forwardModelType)
        if ~isempty(forwardModelType)
            validateattributes(forwardModelType,...
                               {'char'},{'scalartext','nonempty'},...
                               mfilename,'transform',1);
            if ~ismember(forwardModelType,obj.listOfSupportedforwardModelType)
                error(strcat('Variable forwardModelType contains a method that is not supported.\n',...
                             'Supported forward models are: %s.'),...
                             strjoin(cellfun(@(x) sprintf('''%s''',x),obj.listOfSupportedforwardModelType,'UniformOutput',false),', '));
            end
        else
            forwardModelType = 'Identity';
        end
        obj.forwardModelType = forwardModelType;
    end
   function set.noiseType(obj,noiseType)
        if ~isempty(noiseType)
            validateattributes(noiseType,...
                               {'char'},{'scalartext','nonempty'},...
                               mfilename,'transform',1);
            if ~ismember(noiseType,obj.listOfSupportednoiseType)
                error(strcat('Variable noiseType contains a method that is not supported.\n',...
                             'Supported noise types are: %s.'),...
                             strjoin(cellfun(@(x) sprintf('''%s''',x),obj.listOfSupportednoiseType,'UniformOutput',false),', '));
            end
        else
            noiseType = 'Gaussian';
        end
        obj.noiseType = noiseType;
    end

    
end

end