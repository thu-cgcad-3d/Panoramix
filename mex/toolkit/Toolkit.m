
classdef Toolkit < handle
    % properties (SetAccess = private, Hidden = true)
    %     objectHandle; % Handle to the underlying C++ class instance
    % end
    methods (Static)
        function varargout = saveMatToPanoramix(varargin)
             [varargout{1:nargout}] = mexToolkit('saveMatToPanoramix', varargin{:});
        end
		function varargout = loadMatFromPanoramix(varargin)
             [varargout{1:nargout}] = mexToolkit('loadMatFromPanoramix', varargin{:});
        end
    end
end