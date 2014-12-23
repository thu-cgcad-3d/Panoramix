#ifndef PANORAMIX_VIS_VISUALIZERS_HPP
#define PANORAMIX_VIS_VISUALIZERS_HPP

#include "../core/feature.hpp"
#include "../core/cameras.hpp"
#include "../core/misc.hpp"
#include "basic_types.hpp"

namespace panoramix {
    namespace vis {

        struct Options {
            std::string winName;
            Color backgroundColor;
            RenderModeFlags renderMode;
            
            core::PerspectiveCamera camera;

            core::Enabled<ColorTable> colorTable;
            core::Enabled<float> pointSize;
            core::Enabled<float> lineWidth;
            core::Enabled<core::Point3> projectionCenter;
        };


        struct Resource {
            virtual bool isNull() const = 0;
            virtual bool initialize() { return true; }
            virtual bool bind() { return true; }
            virtual ~Resource() {}
        };
        using ResourcePtr = std::unique_ptr<Resource>;


        ResourcePtr MakeTextureResource(const core::Image & tex);


        class VisualObject {
        public:
            VisualObject(const OpenGLShaderSource & shaderSource, VisualObject * parent = nullptr);
            virtual ~VisualObject();

            // initialize rendering
            void initialize() const;

            // render with a given model+view+projection matrix
            void render(const Options & options, const core::Mat4 & mat) const;

            // render with given camera and the stored model matrix
            void renderWithCamera(const Options & options, const core::PerspectiveCamera & cam) const;

            // model matrix
            core::Mat4 & modelMatrix() { return _modelMat; }
            const core::Mat4 & modelMatrix() const { return _modelMat; }

            // mesh
            TriMesh & mesh() { return _mesh; }
            const TriMesh & mesh() const { return _mesh; }

            // resource
            std::vector<ResourcePtr> & resources() { return _resources; }
            const std::vector<ResourcePtr> & resources() const { return _resources; }            


            // traverse its self and its all children
            template <class ConstCallbackFunctorT>
            int depthFirstTraverse(const ConstCallbackFunctorT & callback) const {
                if (!callback(this))
                    return 0;
                int visited = 1;
                for (const VisualObject * ch : _children) {
                    int newlyVisited = ch->depthFirstTraverse(callback);
                    if (newlyVisited == 0)
                        return 0;
                    visited += newlyVisited;
                }
                return visited;
            }

            template <class CallbackFunctorT>
            int depthFirstTraverse(const CallbackFunctorT & callback) {
                if (!callback(this))
                    return 0;
                int visited = 1;
                for (VisualObject * ch : _children) {
                    int newlyVisited = ch->depthFirstTraverse(callback);
                    if (newlyVisited == 0)
                        return 0;
                    visited += newlyVisited;
                }
                return visited;
            }

            // structure
            inline VisualObject * parent() const { return _parent; }
            inline const std::vector<VisualObject*> & children() const { return _children; }

        protected:
            core::Mat4 _modelMat;
            TriMesh _mesh;
            std::vector<ResourcePtr> _resources;
            OpenGLShaderSource _shaderSource;

            VisualObject * _parent;
            std::vector<VisualObject*> _children;

            void * _internal;
        };






        class Visualizer {



        };
        

    }
}
 
#endif