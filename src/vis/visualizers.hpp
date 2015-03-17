#ifndef PANORAMIX_VIS_VISUALIZERS_HPP
#define PANORAMIX_VIS_VISUALIZERS_HPP

#include "../core/feature.hpp"
#include "../core/cameras.hpp"
#include "../core/iterators.hpp"
#include "../core/containers.hpp"
#include "../core/generic_topo.hpp"
#include "basic_types.hpp"

class QWidget;
class QAction;

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
        inline void Discretize(TriMesh & mesh, const core::Point<T, 3> & p, const DiscretizeOptions & o){
            TriMesh::Vertex v;
            v.position = core::Vec4f(p[0], p[1], p[2], 1.0f);
            v.color = o.color;
            v.entityIndex = o.index;
            mesh.addVertex(v);
        }

        template <class T>
        inline void Discretize(TriMesh & mesh, const core::Line<T, 3> & l, const DiscretizeOptions & o){
            TriMesh::Vertex v1, v2;
            v1.position = core::Concat(core::vec_cast<float>(l.first), 1.0f);
            v1.color = o.color;
            v1.entityIndex = o.index;
            v2.position = core::Concat(core::vec_cast<float>(l.second), 1.0f);
            v2.color = o.color;
            v2.entityIndex = o.index;
            mesh.addIsolatedLine(v1, v2);
        }

        void Discretize(TriMesh & mesh, const core::Sphere3 & s, const DiscretizeOptions & o);

        template <class T>
        inline void Discretize(TriMesh & mesh, const core::Sphere<T, 3> & s, const DiscretizeOptions & o){
            Discretize(mesh, core::Sphere3{ vec_cast<double>(s.center), static_cast<double>(s.radius) });
        }

        template <class T>
        inline void Discretize(TriMesh & mesh, const core::Polygon<T, 3> & p, const DiscretizeOptions & o) {
            std::vector<TriMesh::VertHandle> vhandles(p.corners.size());
            for (int i = 0; i < p.corners.size(); i++){
                TriMesh::Vertex v;
                v.position = core::Vec4f(p.corners[i][0], p.corners[i][1], p.corners[i][2], 1.0);
                v.normal = core::vec_cast<float>(p.normal);
                v.color = o.color;
                v.entityIndex = o.index;
                vhandles[i] = mesh.addVertex(v);
            }
            mesh.addPolygon(vhandles);
        }

        void Discretize(TriMesh & mesh, const SpatialProjectedPolygon & spp, const DiscretizeOptions & o);
        
        template <class T, class AllocT>
        inline void Discretize(TriMesh & mesh, const std::vector<T, AllocT> & v, const DiscretizeOptions & o){
            auto oo = o;
            for (auto & e : v){
                Discretize(mesh, e, oo);
                oo.index++;
            }
        }

        template <class T>
        inline void Discretize(TriMesh & mesh, const core::Classified<T> & c, const DiscretizeOptions & o){
            auto oo = o;
            oo.color = o.colorTable[c.claz];
            Discretize(mesh, c.component, oo);
        }

        template <class T>
        inline void Discretize(TriMesh & mesh, const Colored<T> & c, const DiscretizeOptions & o){
            auto oo = o;
            oo.color = c.color;
            Discretize(mesh, c.component, oo);
        }

        template <class AttachedT, class T>
        inline void Discretize(TriMesh & mesh, const std::pair<AttachedT, T> & p, const DiscretizeOptions & o){
            Discretize(mesh, p.second, o);
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
            virtual void initialize() = 0;
            virtual bool isInitialized() const = 0;
            virtual bool bind() = 0;
            virtual void destroy() = 0;
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
            static void clear();
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



        class Visualizer {
        public:
            explicit Visualizer(const std::string & winName = "panoramix::vis::Visualizer") {
                _activeOH = _tree.addRoot(std::make_shared<VisualObject>());

                renderOptions.winName = winName;
                renderOptions.backgroundColor = vis::ColorTag::White;
                renderOptions.renderMode = vis::RenderModeFlag::All;
                renderOptions.camera = core::PerspectiveCamera(500, 500, 250, { 1.0, 1.0, 1.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 });
                renderOptions.bwColor = 0.3;
                renderOptions.bwTexColor = 0.7;
                renderOptions.showInside = true;

                installingOptions.discretizeOptions.color = vis::ColorTag::Black;
                installingOptions.discretizeOptions.colorTable = vis::PredefinedColorTable(vis::ColorTableDescriptor::AllColors);
                installingOptions.discretizeOptions.isolatedTriangles = false;
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

            inline Visualizer & renderMode(vis::RenderModeFlags flags) { renderOptions.renderMode = flags; return *this; }
            inline Visualizer & camera(const core::PerspectiveCamera & cam) { renderOptions.camera = cam; return *this; }

            enum CameraScalePolicy {
                WatchAtMedianScale,
                WatchAtMeanScale,
                WatchAtMaxScale
            };

            QWidget * createWidget(bool autoSetCamera, QWidget * parent);
            void show(bool doModal = true, bool autoSetCamera = true, CameraScalePolicy csp = WatchAtMedianScale);

        private:
            VisualObjectTree _tree;
            VisualObjectHandle _activeOH;
        };

       





    }

}
 
#endif