
classdef Toolkit < handle
    % properties (SetAccess = private, Hidden = true)
    %     objectHandle; % Handle to the underlying C++ class instance
    % end
    methods (Static)
        function varargout = staticMethod(varargin)
             [varargout{1:nargout}] = mexToolkit('showEdges', varargin{:});
        end
    end
end