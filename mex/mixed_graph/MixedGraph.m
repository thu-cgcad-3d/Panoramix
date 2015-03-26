
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
            setpath;
            mexMixedGraph('delete', this.objectHandle);
        end

        function varargout = print(this, varargin)
            setpath;
            [varargout{1:nargout}] = mexMixedGraph('print', this.objectHandle, varargin{:});
        end

        function varargout = show(this, varargin)
            setpath;
            [varargout{1:nargout}] = mexMixedGraph('show', this.objectHandle, varargin{:});
        end

        function varargout = getDepths(this, varargin)
            setpath;
            [varargout{1:nargout}] = mexMixedGraph('getDepths', this.objectHandle, varargin{:});
        end
    end

    methods (Static)
        function varargout = staticMethod(varargin)
            setpath;
        end
    end
end