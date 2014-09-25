#ifndef PANORAMIX_VIS_RENDERABLE_OBJECT_HPP
#define PANORAMIX_VIS_RENDERABLE_OBJECT_HPP

#include "basic_types.hpp"
#include "../core/feature.hpp"

namespace panoramix {
    namespace vis {

        // renderable object
        class RenderableObject {
        public:
            RenderableObject(RenderableObject * parent = nullptr);
            ~RenderableObject();

            // initialize rendering
            virtual void initialize() const {}

            // render with a given model+view+projection matrix
            virtual void render(RenderModeFlags mode, const core::Mat4 & mat) const {}

            // render with given camera and the stored model matrix
            void renderWithCamera(RenderModeFlags mode, const core::PerspectiveCamera & cam) const;

            // bounding box ignoring model matrix
            virtual core::Box3 primaryBoundingBox() const { return core::Box3(); }

            // intersection test ignoring model matrix
            virtual float distanceTo(const core::InfiniteLine3 & ray) const { return std::numeric_limits<float>::max(); }
            virtual bool intersectsWith(const core::InfiniteLine3 & ray, float thres) const { return distanceTo(ray) <= thres; }

            // set texture
            virtual void setTexture(const core::Image & im) {}

            core::Mat4 & modelMatrix() { return _modelMat; }
            const core::Mat4 & modelMatrix() const { return _modelMat; }

            // traverse its self and its all children
            template <class ConstCallbackFunctorT>
            int depthFirstTraverse(ConstCallbackFunctorT && callback) const {
                if (!callback(this))
                    return 0;
                int visited = 1;
                for (const RenderableObject * ch : _children) {
                    int newlyVisited = ch->depthFirstTraverse(callback);
                    if (newlyVisited == 0)
                        return 0;
                    visited += newlyVisited;
                }
                return visited;
            }

            template <class CallbackFunctorT>
            int depthFirstTraverse(CallbackFunctorT && callback) {
                if (!callback(this))
                    return 0;
                int visited = 1;
                for (RenderableObject * ch : _children) {
                    int newlyVisited = ch->depthFirstTraverse(callback);
                    if (newlyVisited == 0)
                        return 0;
                    visited += newlyVisited;
                }
                return visited;
            }

            // structure
            inline RenderableObject * parent() const { return _parent; }
            inline const std::vector<RenderableObject*> children() const { return _children; }

        protected:
            core::Mat4 _modelMat;
            RenderableObject * _parent;
            std::vector<RenderableObject*> _children;
        };


        // default render state
        struct DefaultRenderState {
            DefaultRenderState();
            vis::Color foregroundColor;
            float lineWidth;
            float pointSize;
            vis::ColorTable colorTable;
        };


        // point
        RenderableObject * MakeRenderable(const core::Point3 & p,
            const DefaultRenderState & state = DefaultRenderState(), RenderableObject * parent = nullptr);
        // points
        RenderableObject * MakeRenderable(const std::vector<core::Point3> & points,
            const DefaultRenderState & state = DefaultRenderState(), RenderableObject * parent = nullptr);

        // line
        RenderableObject * MakeRenderable(const core::Line3 & line,
            const DefaultRenderState & state = DefaultRenderState(), RenderableObject * parent = nullptr);
        // lines
        RenderableObject * MakeRenderable(const std::vector<core::Line3> & lines,
            const DefaultRenderState & state = DefaultRenderState(), RenderableObject * parent = nullptr);

        //inline int MakeRenderable(int a, const DefaultRenderState & state = DefaultRenderState(), RenderableObject * parent = nullptr) { return 0; }


        // classified
        template <class T>
        inline RenderableObject * MakeRenderable(const core::Classified<T> & c,
            const DefaultRenderState & state = DefaultRenderState(), RenderableObject * parent = nullptr) {
            auto paramCopy = state;
            paramCopy.foregroundColor = paramCopy.colorTable[c.claz];
            return MakeRenderable(c.component, paramCopy, parent);
        }



        namespace {
            template <class T>
            struct CanMakeRenderableImpl {
                template <class TT>  
                static auto test(int) -> decltype(
                    vis::MakeRenderable(std::declval<TT>(), std::declval<DefaultRenderState>(), nullptr),
                    std::true_type()
                );
                template <class>  
                static std::false_type test(...); 
                static const bool value = std::is_same<decltype(test<T>(0)), std::true_type>::value;
            };
        }

        template <class T>
        struct CanMakeRenderable : public std::integral_constant<bool, CanMakeRenderableImpl<T>::value> {};


        //// opengl object
        //class OpenGLObject : public RenderableObject {
        //public:
        //    explicit OpenGLObject(RenderableObject * parent = nullptr);
        //    ~ OpenGLObject();

        //    void setUpShaders(const OpenGLShaderSource & ss);
        //    void setUpMesh(const OpenGLMesh & mesh);
        //    void setUpTexture(const QImage & tex);

        //    virtual void render(RenderModeFlags mode, const QMatrix4x4 & mat) const override;

        //protected:
        //    void error(const QString & message);

        //private:
        //    OpenGLMesh _mesh;
        //    QOpenGLShaderProgram * _program;
        //    QOpenGLTexture * _texture;
        //};






    }

}
 
#endif