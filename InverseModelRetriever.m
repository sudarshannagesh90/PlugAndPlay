classdef InverseModelRetriever < matlab.mixin.Copyable
% CLASS INVERSEMODELGENERATOR -  Generates the inverse-model coherent
%                                images under the assumption of the given 
%                                forward-model. 
%                                  .
%              
% Constructor:
%   obj = InverseModelRetriever(objectSizePixels)
%
% Inputs:
%    objectSizePixels          : Size of the object in pixels given as [rows,cols]. 
%    inversionModelType        : Type of inversion-model: Can be ML, EM, EMPnP.
%    noiseType                 : Noise type of the measurement (can be one of 'Poisson' or 'Gaussian').
%    sigma_w                   : Std. deviation of the additive
%    white-Gaussian noise (optional). 
%
%
% Output:
%    obj                       : InverseModelRetriever type object with 
%                              overloaded functions.
%
% Overloaded Methods:
%    conj(obj)      : Returns a copy of the object.
%    transpose(obj) : Returns a copy of the operator.
%    ctranspose(obj): Returns a copy of the operator.
%    mtimes(obj,x)  : Called when used as (obj*x). Applies the specified 
%                     InverseModelRetriever on the measurements x after reshaping it to
%                     objectSizePixels. The number of elements in x must be
%                     prod(objectSizePixels).
%

% Author   : Sudarshan Nagesh             
% Institute: NorthWestern University (NU) 
properties(Constant)
    listOfSupportedinversionModelType  = {'ML','EM','EMPnP'};
    listOfSupportednoiseType        = {'Poisson','Gaussian'};
end


properties
    % Properties provided as input.
    objectSizePixels
    inversionModelType
    noiseType
    sigma_w
end

methods
    % Function detectorSamplingOperator - Constructor.
    function obj = InverseModelRetriever(objectSizePixels,inversionModelType,noiseType,sigma_w)
        obj.objectSizePixels                       = objectSizePixels;
        obj.inversionModelType                     = inversionModelType;
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
        x   = reshape(x,obj.objectSizePixels);
        if (strcmp(obj.noiseType,'Gaussian'))
            if (strcmp(obj.inversionModelType,'ML'))
                r        = abs(x).^2-obj.sigma_w^2;
            end
        end
        res = [r(:)];
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
    function set.inversionModelType(obj,inversionModelType)
        if ~isempty(inversionModelType)
            validateattributes(inversionModelType,...
                               {'char'},{'scalartext','nonempty'},...
                               mfilename,'transform',1);
            if ~ismember(inversionModelType,obj.listOfSupportedinversionModelType)
                error(strcat('Variable inversionModelType contains a method that is not supported.\n',...
                             'Supported inversion model types are: %s.'),...
                             strjoin(cellfun(@(x) sprintf('''%s''',x),obj.listOfSupportedinversionModelType,'UniformOutput',false),', '));
            end
        else
            inversionModelType = 'ML';
        end
        obj.inversionModelType = inversionModelType;
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