%--------------------------------------------------------------------------
% General description of an electric field pulse as Matlab objects
%--------------------------------------------------------------------------

classdef generic < handle
    
    properties (Access = public)
        
        fwhm        % FWHM: full width at half maximum
        length      % Full duration of pulse
        delay       % Time delay (center of pulse)
        
        ampli       % Field amplitude (scalar or row vector)
        
        frequ       % Frequency (photon energy)
        linear      % Linear chirp
        quadratic   % Quadratic chirp
        phase       % Phase shift
        
        ind         % index of this pulse
        
    end
    
    properties  (Access = private)
        l_cycle     % length of optical cycles
        n_cycle     % number of optical cycles
        
    end
    
    methods  (Access = public)
        
        % Constructor: Set default values
        function obj = generic
            
            obj.fwhm = 0;
            obj.length = 0;
            obj.delay = 0;
            obj.frequ = 0;
            obj.linear = 0;
            obj.quadratic = 0;
            obj.phase = 0;
        end
        
        % Initialize/check some properties of the pulse
        function init ( obj )
            if obj.frequ~=0
                obj.l_cycle = 2*pi/obj.frequ;
                obj.n_cycle = obj.length / obj.l_cycle;
            end
            if ~isrow(obj.ampli)
                prt.error ('Amplitude of pulse should be given as a scalar or row vector')
            end
        end
        
        % Display pulse properties, thus overloading default 'disp' method
        function disp(obj)
            prt.disp ( ' ' )
            prt.disp (['FWHM duration of pulse      : ' num2str(obj.fwhm)])
            prt.disp (['Full duration of pulse      : ' num2str(obj.length)])
            prt.disp (['Time delay (center of pulse): ' num2str(obj.delay)])
            prt.disp ( ' ' )
            prt.disp (['Field amplitude (vector)    : ' num2str(obj.ampli)])
            prt.disp ( ' ' )
            prt.disp (['Frequency (photon energy)   : ' num2str(obj.frequ)])
            prt.disp (['Linear chirp                : ' num2str(obj.linear)])
            prt.disp (['Quadratic chirp             : ' num2str(obj.quadratic)])
            prt.disp (['Phase shift                 : ' num2str(obj.phase)])
            if obj.frequ~=0
                prt.disp ( ' ' )
                prt.disp (['Duration of optical cycle   : ' num2str(obj.l_cycle)])
                prt.disp (['Number of optical cycles    : ' num2str(obj.n_cycle)])
            end
        end
        
        % Calculate the oscillatory portion of the electric field
        % Constant or (linearly or quadratically) chirped frequencies
        function oscillate = oscillate (obj,timesteps)
            tdelayed = timesteps - obj.delay;
            omega = obj.frequ * ones(size(tdelayed));
            omega = omega + obj.linear * tdelayed + obj.quadratic/2 * tdelayed.^2;
            oscillate = cos ( omega .* tdelayed + obj.phase );
        end
        
    end
    
end


