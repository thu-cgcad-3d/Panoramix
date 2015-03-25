
classdef BoundaryGraph < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance 
        function this = BoundaryGraph(varargin)
            this.objectHandle = mexBoundaryGraph('new', varargin{:});
        end
        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            mexBoundaryGraph('delete', this.objectHandle);
        end

        function varargout = print(this, varargin)
            [varargout{1:nargout}] = mexBoundaryGraph('print', this.objectHandle, varargin{:});
        end

    end
end