#ifndef PANORAMIX_VIS_VISUALIZERS_HPP
#define PANORAMIX_VIS_VISUALIZERS_HPP

#include "../core/feature.hpp"
#include "../core/cameras.hpp"
#include "../core/misc.hpp"
#include "../core/containers.hpp"
#include "../core/graph.hpp"
#include "basic_types.hpp"

namespace panoramix {
    namespace vis {

        enum InteractionID {
            ClickLeftButton,
            PressSpace
        };

        struct Options {
            std::string winName;

            Color backgroundColor;
            RenderModeFlags renderMode;
            
            core::PerspectiveCamera camera;

            float bwColor;
            float bwTexColor;

            float pointSize;
            float lineWidth;
        };


        struct Resource {
            virtual bool isNull() const = 0;
            virtual bool initialize() { return true; }
            virtual bool bind() { return true; }
            virtual ~Resource() {}
        };
        
        using ResourcePtr = std::shared_ptr<Resource>;
        ResourcePtr MakeResource(const core::Image & tex);

        struct ResourceStore {
            static void set(const std::string & name, ResourcePtr r);
            template <class T>
            static void set(const std::string & name, const T & data) { set(name, MakeResource(data)); }
            static ResourcePtr get(const std::string & name);
            static bool has(const std::string & name);
        };

       

        class VisualObject;
        using VisualObjectPtr = std::shared_ptr<VisualObject>;
        using VisualObjectTree = core::Forest<VisualObjectPtr>;

        using VisualObjectHandle = VisualObjectTree::NodeHandle;
        using VisualObjectMeshTriangle = std::pair<VisualObjectHandle, TriMesh::TriangleHandle>;

        typedef bool VisualObjectCallbackFunction(InteractionID, const VisualObjectTree &, const VisualObjectMeshTriangle &);
        
        template <class T>
        using VisualObjectCallbackFunctionSimple = void(InteractionID, T & data);

        

        class VisualObjectInternal;
        class VisualObject {
        public:
            explicit VisualObject();
            explicit VisualObject(const OpenGLShaderSource & shaderSource);
            virtual ~VisualObject();

            // install shaders
            // before initialize()
            void setShaderSource(const OpenGLShaderSource & shaderSource);

            // initialize rendering
            void initialize() const;

            // render with given camera 
            void render(const Options & options, const core::Mat4 & thisModelMatrix) const;

            // model matrix
            core::Mat4 & modelMatrix() { return _modelMat; }
            const core::Mat4 & modelMatrix() const { return _modelMat; }

            // mesh
            TriMesh & mesh() { return _mesh; }
            const TriMesh & mesh() const { return _mesh; }

            // resource
            std::vector<ResourcePtr> & resources() { return _resources; }
            const std::vector<ResourcePtr> & resources() const { return _resources; }

            // data binding   
            template <class FunT>
            inline void bindCallbackFunction(const FunT & fun) { _callback = fun; }
            template <class FunT>
            inline void bindCallbackFunction(FunT && fun) { _callback = std::move(fun); }

            inline bool invokeCallbackFunction(InteractionID iid,
                const VisualObjectTree & tree, 
                const VisualObjectMeshTriangle & omt) const {
                if (_callback)
                    return _callback(iid, tree, omt); 
                return false;
            }

        public:
            struct Flags {
                bool isSelected;
            } flags;

        protected:
            core::Mat4 _modelMat;
            TriMesh _mesh;
            core::Point3 _projectionCenter;
            std::vector<ResourcePtr> _resources;
            OpenGLShaderSource _shaderSource;

            VisualObjectInternal * _internal;
            
            std::function<VisualObjectCallbackFunction> _callback;
        };


        template <class T>
        inline std::shared_ptr<VisualObject> Visualize(const T & data,
            const DiscretizeOptions & dopt = DiscretizeOptions()){
            std::shared_ptr<VisualObject> vo = std::make_shared<VisualObject>();
            Discretize(vo->mesh(), data, dopt);
            vo->setShaderSource(dopt.defaultShaderSource);
            return vo;
        }

        template <class T, class FunT>
        inline std::shared_ptr<VisualObject> VisualizeWithBindingDetail(T & data, const FunT & fun,
            const DiscretizeOptions & dopt = DiscretizeOptions()){
            std::shared_ptr<VisualObject> vo = std::make_shared<VisualObject>();
            Discretize(vo->mesh(), data, dopt);
            vo->setShaderSource(dopt.defaultShaderSource);
            vo->bindCallbackFunction(fun);
            return vo;
        }

        template <class T, class FunT>
        inline std::shared_ptr<VisualObject> VisualizeWithBinding(T & data, const FunT & fun,
            const DiscretizeOptions & dopt = DiscretizeOptions()){
            std::shared_ptr<VisualObject> vo = std::make_shared<VisualObject>();
            Discretize(vo->mesh(), data, dopt);
            vo->setShaderSource(dopt.defaultShaderSource);

            struct CallbackFunction {
                inline CallbackFunction(T & d, const FunT & f)
                    : originalData(d), originalFun(f){}
                inline bool operator() (InteractionID iid, 
                    const VisualObjectTree &, const VisualObjectMeshTriangle &) const {
                    originalFun(iid, originalData);
                    return true;
                }
                T & originalData;
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
                clear();
                _tree = tree;
                update();
            }

            inline void install(VisualObjectTree && tree) {
                clear();
                _tree = std::move(tree);
                update();
            }

            inline const VisualObjectTree & tree() const { return _tree; }

            inline void select(VisualObjectHandle oh) { _tree.data(oh)->flags.isSelected = true; }
            inline void clearSelection() {
                for (auto & n : _tree.nodes())
                    n.data->flags.isSelected = false;
            }

            void update();
            void clear();

            core::Box3 boundingBox() const;
            core::Box3 boundingBoxOfObject(VisualObject * ro) const;
            core::Box3 boundingBoxOfTriangleInObjectMesh(const VisualObjectMeshTriangle & omt) const;

            void initialize() const;
            void render(const Options & options) const;

            VisualObjectMeshTriangle pickOnScreen(const Options & options,
                const core::Point2 & pOnScreen, double distThresOnScreen = 1.0) const;

        private:
            VisualObjectSceneInternal * _internal;
            VisualObjectTree _tree;
        };


        class Visualizer {
        public:
            explicit Visualizer(const std::string & winName = "panoramix::vis::Visualizer") {
                _activeOH = _tree.addRoot(std::make_shared<VisualObject>());

                options.winName = winName;
                options.backgroundColor = ColorFromTag(vis::ColorTag::White);
                options.renderMode = vis::RenderModeFlag::All;
                options.camera = core::PerspectiveCamera(500, 500, 250, { 1.0, 1.0, 1.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 });
                options.bwColor = 0.3;
                options.bwTexColor = 0.7;
                options.pointSize = 10.0;
                options.lineWidth = 5.0;

                doptions.color = vis::ColorFromTag(vis::ColorTag::Black);
                doptions.colorTable = vis::PredefinedColorTable(vis::ColorTableDescriptor::AllColors);
                doptions.isolatedTriangles = true;
                doptions.subdivisionNums[0] = 32;
                doptions.subdivisionNums[1] = 64;
                doptions.defaultShaderSource = vis::PredefinedShaderSource(vis::OpenGLShaderSourceDescriptor::XTriangles);
            }

            Options options;
            DiscretizeOptions doptions;

        public:
            const VisualObjectTree & tree() const { return _tree; }

            VisualObjectHandle activeObjectHandle() const { return _activeOH; }
            VisualObject & activeObject() const { return *_tree.data(_activeOH); }

            const TriMesh & activeMesh() const { return _tree.data(_activeOH)->mesh(); }
            TriMesh & activeMesh() { return _tree.data(_activeOH)->mesh(); }

        public:

            template <class T>
            inline Visualizer & add(const T & data) {
                _tree.add(_activeOH, Visualize<T>(data, doptions));
                return *this;
            }
            template <class T, class FunT>
            inline Visualizer & addWithBindingDetail(T & data, const FunT & fun) {
                _tree.add(_activeOH, VisualizeWithBindingDetail<T, FunT>(data, fun, doptions));
                return *this;
            }
            template <class T, class FunT>
            inline Visualizer & addWithBinding(T & data, const FunT & fun) {
                _tree.add(_activeOH, VisualizeWithBinding<T, FunT>(data, fun, doptions));
                return *this;
            }


            template <class T>
            inline Visualizer & begin(const T & data) {
                _activeOH = _tree.add(_activeOH, Visualize<T>(data, doptions));
                return *this;
            }
            template <class T, class FunT>
            inline Visualizer & beginWithBindingDetail(T & data, const FunT &fun) {
                _activeOH = _tree.add(_activeOH, VisualizeWithBindingDetail<T, FunT>(data, fun, doptions));
                return *this;
            }
            template <class T, class FunT>
            inline Visualizer & beginWithBinding(T & data, const FunT & fun) { 
                _activeOH = _tree.add(_activeOH, VisualizeWithBinding<T, FunT>(data, fun, doptions));
                return *this;
            }


            inline Visualizer & shaderSource(const OpenGLShaderSource & ss) { 
                activeObject().setShaderSource(ss);
                return *this; 
            }
            inline Visualizer & resource(const std::string resourceName) {
                activeObject().resources().push_back(ResourceStore::get(resourceName));
                return *this;
            }


            inline Visualizer & end() { 
                if (!_tree.isRoot(_activeOH)){
                    _activeOH = _tree.parent(_activeOH);
                }
                return *this; 
            }


            void show(bool doModal = true, bool autoSetCamera = true);

        private:
            VisualObjectTree _tree;
            VisualObjectHandle _activeOH;
        };

       





    }
}
 
#endif