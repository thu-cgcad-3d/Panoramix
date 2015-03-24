#ifndef PANORAMIX_SCENE_HPP
#define PANORAMIX_SCENE_HPP

#include "../core/generic_topo.hpp"
#include "../core/utilities.hpp"
#include "../core/cameras.hpp"
#include "basic_types.hpp"
#include "discretization.hpp"
#include "resource.hpp"
 
namespace panoramix {
    namespace gui {


        enum InteractionID {
            ClickLeftButton,
            PressSpace,
            Unknown
        };


 
        // building visual objects
        class VisualObject;
        using VisualObjectPtr = std::shared_ptr<VisualObject>;
        using VisualObjectTree = core::Forest<VisualObjectPtr>;

        using VisualObjectHandle = VisualObjectTree::NodeHandle;
        using VisualObjectMeshTriangle = std::pair<VisualObjectHandle, TriMesh::TriangleHandle>;
        using VisualObjectEntityID = std::pair<VisualObjectHandle, int>;

        typedef bool VisualObjectCallbackFunction(InteractionID, const VisualObjectTree &, const VisualObjectEntityID &);
        template <class T> using VisualObjectCallbackFunctionSimple = void(InteractionID, T & data);

        namespace {

            template <class FunAndArgT>
            struct IsComplexCallbackFunctionImp {
                template <class FunAndArgTT>
                static auto test(int) -> decltype(
                    std::declval<FunAndArgTT>().first(Unknown,
                    std::declval<VisualObjectTree>(),
                    std::declval<VisualObjectEntityID>()),
                    std::true_type()
                    );

                template <class>
                static std::false_type
                    test(...);

                enum { value = std::is_same<decltype(test<FunAndArgT>(0)), std::true_type>::value };
            };

            template <class FunAndArgT>
            struct IsSimpleCallbackFunctionImp {
                template <class FunAndArgTT>
                static auto test(int) -> decltype(
                    std::declval<FunAndArgTT &>().first(Unknown, std::declval<FunAndArgTT &>().second),
                    std::true_type()
                    );

                template <class>
                static std::false_type
                    test(...);

                enum { value = std::is_same<decltype(test<FunAndArgT>(0)), std::true_type>::value };
            };

        }


        enum class CallbackFunctionType {
            Complex, Simple, None
        };
        using ComplexTag = std::integral_constant<CallbackFunctionType, CallbackFunctionType::Complex>;
        using SimpleTag = std::integral_constant<CallbackFunctionType, CallbackFunctionType::Simple>;
        using NoneTag = std::integral_constant<CallbackFunctionType, CallbackFunctionType::None>;

        template <class FunT, class T>
        struct CallbackFunctionTraits
            : std::integral_constant<CallbackFunctionType,
            (IsComplexCallbackFunctionImp<std::pair<FunT, T>>::value ? CallbackFunctionType::Complex :
            (IsSimpleCallbackFunctionImp<std::pair<FunT, T>>::value ? CallbackFunctionType::Simple :
            CallbackFunctionType::None))> {};




        struct RenderOptions {
            std::string winName;
            Color backgroundColor;
            RenderModeFlags renderMode;
            core::PerspectiveCamera camera;
            float bwColor;
            float bwTexColor;
            bool showInside;
        };

        class VisualObjectInternal;
        class VisualObject {
        public:
            explicit VisualObject(int eleNum = 1);
            explicit VisualObject(int eleNum, const OpenGLShaderSource & shaderSource);
            virtual ~VisualObject();

            // install shaders
            // before initialize()
            void setShaderSource(const OpenGLShaderSource & shaderSource);

            // initialize rendering
            void initialize() const;

            // render with given camera 
            void render(const RenderOptions & options, const core::Mat4f & thisModelMatrix) const;

            // model matrix
            core::Mat4f & modelMatrix() { return _modelMat; }
            const core::Mat4f & modelMatrix() const { return _modelMat; }

            // mesh
            TriMesh & mesh() { return _mesh; }
            const TriMesh & mesh() const { return _mesh; }

            bool isSingleEntityMeshTriangle(TriMesh::TriangleHandle t) const;
            int entityIDOfMeshTriangle(TriMesh::TriangleHandle t) const;

            void selectEntity(int entId);
            void switchEntitySelection(int entId);
            void clearSelection();
            inline const std::set<int> & selectedEntities() const { return _selectedEntities; }
            inline bool entityIsSelected(int entId) const { return core::Contains(_selectedEntities, entId); }

            // resource
            std::vector<ResourcePtr> & resources() { return _resources; }
            const std::vector<ResourcePtr> & resources() const { return _resources; }

            // data binding   
            template <class FunT>
            inline void bindCallbackFunction(FunT && fun) { _callback = std::forward<FunT>(fun); }

            inline bool invokeCallbackFunction(InteractionID iid,
                const VisualObjectTree & tree,
                const VisualObjectEntityID & entity) const {
                if (_callback)
                    return _callback(iid, tree, entity);
                return false;
            }

            inline void setLineWidth(float lw){ _lineWidth = lw; }
            inline void setPointSize(float ps){ _pointSize = ps; }

        protected:
            core::Mat4f _modelMat;
            TriMesh _mesh;
            core::Point3 _projectionCenter;
            std::vector<ResourcePtr> _resources;
            OpenGLShaderSource _shaderSource;

            float _lineWidth;
            float _pointSize;

            VisualObjectInternal * _internal;

            std::function<VisualObjectCallbackFunction> _callback;

            int _entitiNum;
            std::set<int> _selectedEntities;
        };


        struct VisualObjectInstallingOptions {
            OpenGLShaderSource defaultShaderSource;
            float lineWidth;
            float pointSize;
            DiscretizeOptions discretizeOptions;
        };


        template <class T>
        std::shared_ptr<VisualObject> Visualize(const T & data,
            const VisualObjectInstallingOptions & o){
            std::shared_ptr<VisualObject> vo = std::make_shared<VisualObject>(1);

            auto dopt = o.discretizeOptions;
            dopt.index = 0;
            Discretize(vo->mesh(), data, dopt);

            vo->setShaderSource(o.defaultShaderSource);
            vo->setLineWidth(o.lineWidth);
            vo->setPointSize(o.pointSize);
            return vo;
        }

        template <class T>
        std::shared_ptr<VisualObject> VisualizeCollection(const std::vector<T> & data,
            const VisualObjectInstallingOptions & o){
            std::shared_ptr<VisualObject> vo = std::make_shared<VisualObject>(data.size());

            auto dopt = o.discretizeOptions;
            dopt.index = 0;
            Discretize(vo->mesh(), data, dopt);

            vo->setShaderSource(o.defaultShaderSource);
            vo->setLineWidth(o.lineWidth);
            vo->setPointSize(o.pointSize);
            return vo;
        }

        template <class T, class FunT>
        std::shared_ptr<VisualObject> Visualize(T & data, const FunT & fun,
            const VisualObjectInstallingOptions & o,
            ComplexTag){

            std::shared_ptr<VisualObject> vo = std::make_shared<VisualObject>(1);

            auto dopt = o.discretizeOptions;
            dopt.index = 0;
            Discretize(vo->mesh(), data, dopt);

            vo->setShaderSource(o.defaultShaderSource);
            vo->bindCallbackFunction(fun);
            vo->setLineWidth(o.lineWidth);
            vo->setPointSize(o.pointSize);
            return vo;
        }

        template <class T, class FunT>
        std::shared_ptr<VisualObject> VisualizeCollection(std::vector<T> & data, const FunT & fun,
            const VisualObjectInstallingOptions & o,
            ComplexTag){

            std::shared_ptr<VisualObject> vo = std::make_shared<VisualObject>(data.size());

            auto dopt = o.discretizeOptions;
            dopt.index = 0;
            Discretize(vo->mesh(), data, dopt);

            vo->setShaderSource(o.defaultShaderSource);
            vo->bindCallbackFunction(fun);
            vo->setLineWidth(o.lineWidth);
            vo->setPointSize(o.pointSize);
            return vo;
        }


        template <class T, class FunT>
        std::shared_ptr<VisualObject> Visualize(T & data, const FunT & fun,
            const VisualObjectInstallingOptions & o,
            SimpleTag){

            std::shared_ptr<VisualObject> vo = std::make_shared<VisualObject>(1);

            auto dopt = o.discretizeOptions;
            dopt.index = 0;
            Discretize(vo->mesh(), data, dopt);

            vo->setShaderSource(o.defaultShaderSource);
            vo->setLineWidth(o.lineWidth);
            vo->setPointSize(o.pointSize);

            struct CallbackFunction {
                inline CallbackFunction(T & d, const FunT & f)
                : originalData(d), originalFun(f){}
                inline bool operator() (InteractionID iid,
                    const VisualObjectTree &, const VisualObjectEntityID &) const {
                    originalFun(iid, originalData);
                    return true;
                }
                T & originalData;
                const FunT & originalFun;
            };

            vo->bindCallbackFunction(CallbackFunction(data, fun));
            return vo;
        }

        template <class T, class FunT>
        std::shared_ptr<VisualObject> VisualizeCollection(std::vector<T> & data, const FunT & fun,
            const VisualObjectInstallingOptions & o,
            SimpleTag){

            std::shared_ptr<VisualObject> vo = std::make_shared<VisualObject>(data.size());

            auto dopt = o.discretizeOptions;
            dopt.index = 0;
            Discretize(vo->mesh(), data, dopt);

            vo->setShaderSource(o.defaultShaderSource);
            vo->setLineWidth(o.lineWidth);
            vo->setPointSize(o.pointSize);

            struct CallbackFunction {
                inline CallbackFunction(std::vector<T> & d, const FunT & f)
                : originalData(d), originalFun(f){}
                inline bool operator() (InteractionID iid,
                    const VisualObjectTree & tree, const VisualObjectEntityID & entityID) const {
                    originalFun(iid, originalData[entityID.second]);
                    return true;
                }
                std::vector<T> & originalData;
                const FunT & originalFun;
            };

            vo->bindCallbackFunction(CallbackFunction(data, fun));
            return vo;
        }






        class VisualObjectSceneInternal;
        class VisualObjectScene {
        public:
            VisualObjectScene();
            explicit VisualObjectScene(const VisualObjectTree & tree);
            ~VisualObjectScene();

            inline void install(const VisualObjectTree & tree) {
                _tree = tree;
                update();
            }

            inline void install(VisualObjectTree && tree) {
                _tree = std::move(tree);
                update();
            }

            inline const VisualObjectTree & tree() const { return _tree; }

            inline void select(VisualObjectEntityID ent) { _tree.data(ent.first)->selectEntity(ent.second); }
            inline void switchSelect(VisualObjectEntityID ent) { _tree.data(ent.first)->switchEntitySelection(ent.second); }
            inline void clearSelection() {
                for (auto & n : _tree.nodes())
                    n.data->clearSelection();
            }

            void update();
            void clear();

            const core::Box3 & boundingBox() const;
            core::Box3 boundingBoxOfObject(VisualObjectHandle h) const;
            core::Box3 boundingBoxOfTriangleInObjectMesh(const VisualObjectMeshTriangle & omt) const;

            void initialize() const;
            void render(const RenderOptions & options) const;

            VisualObjectMeshTriangle pickOnScreen(const RenderOptions & options,
                const core::Point2 & pOnScreen) const;

        private:
            VisualObjectSceneInternal * _internal;
            VisualObjectTree _tree;
        };

 
    }
}
 
#endif