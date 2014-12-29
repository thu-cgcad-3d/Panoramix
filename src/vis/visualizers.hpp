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
            PressSpace,
            Unknown
        };


        // discretization
        struct DiscretizeOptions {
            inline DiscretizeOptions() 
                : color(0, 0, 0, 1), index(0), isolatedTriangles(false) {
                subdivisionNums[0] = 32;
                subdivisionNums[1] = 64;
            }
            Color color;
            vis::ColorTable colorTable;
            int index;
            bool isolatedTriangles;
            int subdivisionNums[2];
        };

        template <class T>
        inline TriMesh & Discretize(TriMesh & mesh, const core::Point<T, 3> & p, const DiscretizeOptions & o){
            auto vh = mesh.addVertex(core::Vec4f(p[0], p[1], p[2], 1.0f));
            mesh.vertices[vh].color = o.color;
            mesh.vertices[vh].entityIndex = o.index;
            return mesh;
        }

        template <class T>
        inline TriMesh & Discretize(TriMesh & mesh, const core::Line<T, 3> & l, const DiscretizeOptions & o){
            TriMesh::Vertex v1, v2;
            v1.position = core::Concat(core::ConvertTo<float>(l.first), 1.0f);
            v1.color = o.color;
            v1.entityIndex = o.index;
            v2.position = core::Concat(core::ConvertTo<float>(l.second), 1.0f);
            v2.color = o.color;
            v2.entityIndex = o.index;
            mesh.addIsolatedLine(v1, v2);
            return mesh;
        }

        TriMesh & Discretize(TriMesh & mesh, const core::Sphere3 & s, const DiscretizeOptions & o);

        template <class T>
        inline TriMesh & Discretize(TriMesh & mesh, const core::Sphere<T, 3> & s, const DiscretizeOptions & o){
            return Discretize(mesh, core::Sphere3{ ConvertTo<double>(s.center), static_cast<double>(s.radius) });
        }

        template <class T>
        inline TriMesh & Discretize(TriMesh & mesh, const core::Polygon<T, 3> & p, const DiscretizeOptions & o) {
            std::vector<TriMesh::VertHandle> vhandles(p.corners.size());
            for (int i = 0; i < p.corners.size(); i++){
                TriMesh::Vertex v;
                v.position = core::Vec4f(p.corners[i][0], p.corners[i][1], p.corners[i][2], 1.0);
                v.normal = core::ConvertTo<float>(p.normal);
                v.color = o.color;
                v.entityIndex = o.index;
                vhandles[i] = mesh.addVertex(v);
            }
            mesh.addPolygon(vhandles);
            return mesh;
        }

        TriMesh & Discretize(TriMesh & mesh, const SpatialProjectedPolygon & spp, const DiscretizeOptions & o);
        
        template <class T, class AllocT>
        inline TriMesh & Discretize(TriMesh & mesh, const std::vector<T, AllocT> & v, const DiscretizeOptions & o){
            auto oo = o;
            for (auto & e : v){
                Discretize(mesh, e, oo);
                oo.index++;
            }
            return mesh;
        }

        template <class T>
        inline TriMesh & Discretize(TriMesh & mesh, const core::Classified<T> & c, const DiscretizeOptions & o){
            auto oo = o;
            oo.color = o.colorTable[c.claz];
            return Discretize(mesh, c.component, oo);
        }


        // Is discretizable ?
        namespace {
            template <class T>
            struct IsDiscretizableImp {
                template <class TT>
                static auto test(int) -> decltype(
                    vis::Discretize(std::declval<TriMesh &>(), std::declval<TT>(), std::declval<DiscretizeOptions>()),
                    std::true_type()
                    );
                template <class>
                static std::false_type test(...);
                static const bool value = std::is_same<decltype(test<T>(0)), std::true_type>::value;
            };
        }

        template <class T>
        struct IsDiscretizable : std::integral_constant<bool, IsDiscretizableImp<T>::value> {};


        // resource making
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

       


        // building visual objects
        class VisualObject;
        using VisualObjectPtr = std::shared_ptr<VisualObject>;
        using VisualObjectTree = core::Forest<VisualObjectPtr>;

        using VisualObjectHandle = VisualObjectTree::NodeHandle;
        using VisualObjectMeshTriangle = std::pair<VisualObjectHandle, TriMesh::TriangleHandle>;

        typedef bool VisualObjectCallbackFunction(InteractionID, 
            const VisualObjectTree &, const VisualObjectMeshTriangle &);        
        template <class T> using VisualObjectCallbackFunctionSimple = void(InteractionID, T & data);

        namespace {

            template <class FunAndArgT>
            struct IsComplexCallbackFunctionImp {
                template <class FunAndArgTT>
                static auto test(int) -> decltype(
                    std::declval<FunAndArgTT>().first(Unknown,
                        std::declval<VisualObjectTree>(), 
                        std::declval<VisualObjectMeshTriangle>()),
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

            inline void setLineWidth(float lw){ _lineWidth = lw; }
            inline void setPointSize(float ps){ _pointSize = ps; }

        public:
            struct Flags {
                bool isSelected;
            } flags;

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
            int _eleNum;
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
            ComplexTag ){
            
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
            ComplexTag ){

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
                    const VisualObjectTree & tree, const VisualObjectMeshTriangle & omt) const {
                    auto & mesh = tree.data(omt.first)->mesh();
                    TriMesh::TriangleHandle v1, v2, v3;
                    mesh.fetchTriangleVerts(omt.second, v1, v2, v3);
                    int indices[] = { 
                        mesh.vertices[v1].entityIndex, 
                        mesh.vertices[v2].entityIndex, 
                        mesh.vertices[v3].entityIndex 
                    };
                    originalFun(iid, originalData[indices[0]]);
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

            inline void select(VisualObjectHandle oh) { _tree.data(oh)->flags.isSelected = true; }
            inline void clearSelection() {
                for (auto & n : _tree.nodes())
                    n.data->flags.isSelected = false;
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


        class Visualizer {
        public:
            explicit Visualizer(const std::string & winName = "panoramix::vis::Visualizer") {
                _activeOH = _tree.addRoot(std::make_shared<VisualObject>());

                renderOptions.winName = winName;
                renderOptions.backgroundColor = ColorFromTag(vis::ColorTag::White);
                renderOptions.renderMode = vis::RenderModeFlag::All;
                renderOptions.camera = core::PerspectiveCamera(500, 500, 250, { 1.0, 1.0, 1.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 });
                renderOptions.bwColor = 0.3;
                renderOptions.bwTexColor = 0.7;

                installingOptions.discretizeOptions.color = vis::ColorFromTag(vis::ColorTag::Black);
                installingOptions.discretizeOptions.colorTable = vis::PredefinedColorTable(vis::ColorTableDescriptor::AllColors);
                installingOptions.discretizeOptions.isolatedTriangles = true;
                installingOptions.discretizeOptions.subdivisionNums[0] = 32;
                installingOptions.discretizeOptions.subdivisionNums[1] = 64;
                installingOptions.defaultShaderSource = vis::PredefinedShaderSource(vis::OpenGLShaderSourceDescriptor::XTriangles);
                installingOptions.pointSize = 10.0;
                installingOptions.lineWidth = 5.0;
            }
                        
            VisualObjectInstallingOptions installingOptions;
            RenderOptions renderOptions;

        public:
            const VisualObjectTree & tree() const { return _tree; }

            VisualObjectHandle activeObjectHandle() const { return _activeOH; }
            VisualObject & activeObject() const { return *_tree.data(_activeOH); }

            const TriMesh & activeMesh() const { return _tree.data(_activeOH)->mesh(); }
            TriMesh & activeMesh() { return _tree.data(_activeOH)->mesh(); }

        public:

            template <class T>
            inline Visualizer & add(const T & data) {
                _tree.add(_activeOH, Visualize(data, installingOptions));
                return *this;
            }
            template <class T, class FunT>
            inline Visualizer & add(T & data, const FunT & fun) {
                _tree.add(_activeOH, Visualize(data, fun, doptions,
                    std::integral_constant<CallbackFunctionType, CallbackFunctionTraits<FunT, T>::value>()));
                return *this;
            }


            template <class T>
            inline Visualizer & begin(const T & data) {
                _activeOH = _tree.add(_activeOH, Visualize<T>(data, installingOptions));
                return *this;
            }

            template <class T>
            inline Visualizer & begin(const std::vector<T> & data) {
                _activeOH = _tree.add(_activeOH, VisualizeCollection<T>(data, installingOptions));
                return *this;
            }

            template <class T, class FunT>
            inline Visualizer & begin(T & data, const FunT & fun) { 
                _activeOH = _tree.add(_activeOH, Visualize<T, FunT>(data, fun, installingOptions,
                    std::integral_constant<CallbackFunctionType, CallbackFunctionTraits<FunT, T>::value>()));
                return *this;
            }

            template <class T, class FunT>
            inline Visualizer & begin(std::vector<T> & data, const FunT & fun) {
                _activeOH = _tree.add(_activeOH, VisualizeCollection<T, FunT>(data, fun, installingOptions,
                    std::integral_constant<CallbackFunctionType, CallbackFunctionTraits<FunT, T>::value>()));
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

            inline Visualizer & camera(const core::PerspectiveCamera & cam) { renderOptions.camera = cam; return *this; }

            void show(bool doModal = true, bool autoSetCamera = true);

        private:
            VisualObjectTree _tree;
            VisualObjectHandle _activeOH;
        };

       





    }

}
 
#endif