
classdef MixedGraph < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance 
        function this = MixedGraph(varargin)
            setpath;
            this.objectHandle = mexMixedGraph('new', varargin{:});
        end
        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            mexMixedGraph('delete', this.objectHandle);
        end

        function varargout = showVPsLines(this, varargin)
            [varargout{1:nargout}] = mexMixedGraph('showVPsLines', this.objectHandle, varargin{:});
        end
        
        function varargout = showROC(this, varargin)
            [varargout{1:nargout}] = mexMixedGraph('showROC', this.objectHandle, varargin{:});
        end

        function varargout = solve(this, varargin)
            [varargout{1:nargout}] = mexMixedGraph('solve', this.objectHandle, varargin{:});
        end

        function varargout = looseC(this, varargin)
            [varargout{1:nargout}] = mexMixedGraph('looseC', this.objectHandle, varargin{:});
        end       

        function varargout = show(this, varargin)
            [varargout{1:nargout}] = mexMixedGraph('show', this.objectHandle, varargin{:});
        end

        function varargout = depths(this, varargin)
            [varargout{1:nargout}] = mexMixedGraph('depths', this.objectHandle, varargin{:});
        end
    end

    methods (Static)
        function varargout = staticMethod(varargin)
        end
    end
end