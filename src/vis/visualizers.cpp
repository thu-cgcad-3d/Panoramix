#include <QtOpenGL>

#include "../core/utilities.hpp"
#include "qt_glue.hpp"
#include "visualizers.hpp"
#include "singleton.hpp"

#define GL_ALPHA_TEST 0x0BC0

namespace panoramix {
    namespace vis {

        using namespace core;




        TriMesh & Discretize(TriMesh & mesh, const core::Sphere3 & s, const DiscretizeOptions & o) {
            int m = o.subdivisionNums[0];
            int n = o.subdivisionNums[1];
            if (!o.isolatedTriangles){
                mesh.vertices.reserve(mesh.vertices.size() + m * n);
                std::vector<std::vector<TriMesh::VertHandle>> vhs(m, std::vector<TriMesh::VertHandle>(n));
                for (int i = 0; i < m; i++) {
                    for (int j = 0; j < n; j++) {
                        float xratio = 1.0f / n * j;
                        float yratio = 1.0f / (m - 1) * i;
                        float xangle = M_PI * 2 * xratio;
                        float yangle = M_PI * yratio - M_PI_2;
                        Vec4 point = {
                            cos(xangle)*cos(yangle) * s.radius + s.center[0],
                            sin(xangle)*cos(yangle) * s.radius + s.center[1],
                            sin(yangle) * s.radius + s.center[2],
                            1
                        };
                        TriMesh::Vertex v;
                        v.position = point;
                        v.texCoord = { xratio, yratio };
                        v.color = Vec4f(o.color[0], o.color[1], o.color[2], 1.0f);
                        v.entityIndex = o.index;
                        vhs[i][j] = mesh.addVertex(v);
                    }
                }
                for (int i = 1; i < m; i++) {
                    int previ = i == 0 ? m - 1 : i - 1;
                    for (int j = 0; j < n; j++) {
                        int prevj = j == 0 ? n - 1 : j - 1;
                        mesh.addTriangle(vhs[i][j], vhs[i][prevj], vhs[previ][prevj]);
                        mesh.addTriangle(vhs[i][j], vhs[previ][prevj], vhs[previ][j]);
                    }
                }
            }
            else {
                std::vector<std::vector<TriMesh::Vertex>> vs(m, std::vector<TriMesh::Vertex>(n));
                for (int i = 0; i < m; i++) {
                    for (int j = 0; j < n; j++) {
                        float xratio = 1.0f / n * j;
                        float yratio = 1.0f / (m - 1) * i;
                        float xangle = M_PI * 2 * xratio;
                        float yangle = M_PI * yratio - M_PI_2;
                        Vec4 point = {
                            cos(xangle)*cos(yangle) * s.radius + s.center[0],
                            sin(xangle)*cos(yangle) * s.radius + s.center[1],
                            sin(yangle) * s.radius + s.center[2],
                            1
                        };
                        TriMesh::Vertex v;
                        v.position = point;
                        v.texCoord = { xratio, yratio };
                        v.color = Vec4f(o.color[0], o.color[1], o.color[2], 1.0f);
                        v.entityIndex = o.index;
                        vs[i][j] = v;
                    }
                }
                for (int i = 1; i < m; i++) {
                    int previ = i == 0 ? m - 1 : i - 1;
                    for (int j = 0; j < n; j++) {
                        int prevj = j == 0 ? n - 1 : j - 1;
                        mesh.addIsolatedTriangle(vs[i][j], vs[i][prevj], vs[previ][prevj]);
                        mesh.addIsolatedTriangle(vs[i][j], vs[previ][prevj], vs[previ][j]);
                    }
                }
            }
            return mesh;
        }



        TriMesh & Discretize(TriMesh & mesh, const SpatialProjectedPolygon & spp, const DiscretizeOptions & o){
            std::vector<Vec3> cs(spp.corners.size());
            for (int i = 0; i < spp.corners.size(); i++){
                InfiniteLine3 line(spp.projectionCenter, spp.corners[i] - spp.projectionCenter);
                cs[i] = IntersectionOfLineAndPlane(line, spp.plane).position;
            }
            std::vector<TriMesh::VertHandle> vhandles(cs.size());
            for (int i = 0; i < cs.size(); i++){
                TriMesh::Vertex v;
                v.position = Concat(cs[i], 1.0);
                v.normal = spp.plane.normal;
                v.color = o.color;
                v.entityIndex = o.index;
                vhandles[i] = mesh.addVertex(v);
            }
            mesh.addPolygon(vhandles);
            return mesh;
        }





        ResourcePtr MakeResource(const core::Image & image){
            
            struct TextureResource : Resource {
                inline TextureResource(const core::Image & im)
                    : initialized(false), 
                    image(MakeQImage(im)), texture(new QOpenGLTexture(QOpenGLTexture::Target2D)) {}
                
                virtual bool isNull() const override { return image.isNull(); }
                virtual void initialize() override {
                    if (!texture->isCreated()){
                        texture->create();
                    }
                    if (!texture->isCreated()){
                        return;
                    }
                    Q_ASSERT(texture->textureId());
                    texture->setData(image.mirrored());
                    texture->setMinificationFilter(QOpenGLTexture::Linear);
                    texture->setMagnificationFilter(QOpenGLTexture::Linear);
                    texture->release();
                    initialized = true;
                }
                virtual bool isInitialized() const override {
                    return initialized;
                }
                virtual void destroy() override {
                    if (texture->isCreated()){
                        texture->destroy();
                    }
                    initialized = false;
                }
                virtual bool bind() override { texture->bind(0); return texture->isBound(); }
                virtual ~TextureResource() { texture->destroy();  delete texture; }
                QImage image;
                QOpenGLTexture * texture;
                bool initialized;
            };

            return std::make_shared<TextureResource>(image);
        }


        static std::unordered_map<std::string, ResourcePtr> g_ResourcesTable;

        void ResourceStore::set(const std::string & name, ResourcePtr r) {
            g_ResourcesTable[name] = r;
        }
        ResourcePtr ResourceStore::get(const std::string & name) {
            if (Contains(g_ResourcesTable, name)){
                return g_ResourcesTable.at(name);
            }
            return nullptr;
        }
        bool ResourceStore::has(const std::string & name){
            return Contains(g_ResourcesTable, name);
        }
        void ResourceStore::clear(){
            g_ResourcesTable.clear();
        }







        namespace {

            template <class T> struct OpenGLDataTraits {};
            template <> struct OpenGLDataTraits<int8_t> { enum { GLType = GL_BYTE }; };
            template <> struct OpenGLDataTraits<uint8_t> { enum { GLType = GL_UNSIGNED_BYTE }; };
            template <> struct OpenGLDataTraits<int16_t> { enum { GLType = GL_SHORT }; };
            template <> struct OpenGLDataTraits<uint16_t> { enum { GLType = GL_UNSIGNED_SHORT }; };
            template <> struct OpenGLDataTraits<int32_t> { enum { GLType = GL_INT }; };
            template <> struct OpenGLDataTraits<uint32_t> { enum { GLType = GL_UNSIGNED_INT }; };
            template <> struct OpenGLDataTraits<float> { enum { GLType = GL_FLOAT }; };
            template <> struct OpenGLDataTraits<double> { enum { GLType = GL_DOUBLE }; };

            template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
            inline void SetAttributeArrayWithTriMeshVertices(
                QOpenGLShaderProgram * program, const char * attributeName,
                const T & firstVector) {
                program->setAttributeArray(attributeName, OpenGLDataTraits<T>::GLType,
                    &firstVector, 1, sizeof(TriMesh::Vertex));
            }

            template <class T, int N, class = std::enable_if_t<std::is_arithmetic<T>::value>>
            inline void SetAttributeArrayWithTriMeshVertices(
                QOpenGLShaderProgram * program, const char * attributeName,
                const core::Vec<T, N> & firstVector) {
                program->setAttributeArray(attributeName, OpenGLDataTraits<T>::GLType,
                    &firstVector, N, sizeof(TriMesh::Vertex));
            }

            template <class T>
            inline void DrawElements(GLenum mode, const std::vector<T> & indices) {
                glDrawElements(mode, indices.size(), OpenGLDataTraits<T>::GLType, indices.data());
            }
        }







        class VisualObjectInternal {
        public:
            QOpenGLShaderProgram * program;            
            inline VisualObjectInternal() :
                program(new QOpenGLShaderProgram) {}
            inline ~VisualObjectInternal(){
                delete program;
            }
        };

        VisualObject::VisualObject(int eleNum)
            : _entitiNum(eleNum),
            _shaderSource(PredefinedShaderSource(OpenGLShaderSourceDescriptor::XLines)),
            _modelMat(core::Mat4::eye()), _projectionCenter(0, 0, 0) {
            _internal = new VisualObjectInternal;
            _lineWidth = 2.0;
            _pointSize = 5.0;
        }

        VisualObject::VisualObject(int eleNum, const OpenGLShaderSource & shaderSource)
            : _entitiNum(eleNum),
            _shaderSource(shaderSource),
            _modelMat(core::Mat4::eye()), _projectionCenter(0, 0, 0) {
            _internal = new VisualObjectInternal;
            _lineWidth = 2.0;
            _pointSize = 5.0;
        }

        VisualObject::~VisualObject(){
            for (auto & res : _resources){
                if (!res->isInitialized())
                    continue;
                res->destroy();
                if (res->isInitialized())
                    qDebug() << (_internal->program->log());
            }

            delete _internal;
        }

        void VisualObject::setShaderSource(const OpenGLShaderSource & shaderSource){
            if (_internal->program->isLinked()){
                qDebug() << "program is already linked! setting shaders failed!";
                return;
            }
            _shaderSource = shaderSource;
        }

        void VisualObject::initialize() const {

            auto vo = _internal;
            auto program = vo->program;

            if (program->isLinked()) {
                return;
            }

            // setup shaders
            if (!program->addShaderFromSourceCode(QOpenGLShader::Vertex,
                QByteArray::fromRawData(_shaderSource.vertexShaderSource().data(), _shaderSource.vertexShaderSource().size())))
                qDebug() << (program->log());
            if (!program->addShaderFromSourceCode(QOpenGLShader::Fragment,
                QByteArray::fromRawData(_shaderSource.fragmentShaderSource().data(), _shaderSource.fragmentShaderSource().size())))
                qDebug() << (program->log());
            if (!program->link() || !program->bind()) {
                qDebug() << (program->log());
                return;
            }
            Q_ASSERT(program->isLinked());

            // initialize resources
            for (auto & res : _resources){
                if (res->isInitialized())
                    continue;
                res->initialize();
                if (!res->isInitialized())
                    qDebug() << (program->log());
            }           

            program->release();

        }

        void VisualObject::render(const RenderOptions & options, const core::Mat4f & thisModelMatrix) const {
            if (_mesh.vertices.empty())
                return;

            auto vo = _internal;
            auto program = vo->program;

            Q_ASSERT(program->isLinked());
            program->bind();

            // bind resources
            for (auto & res : _resources){
                if (!res->bind())
                    qDebug() << (program->log());
            }

            glLineWidth(_lineWidth);
            glPointSize(_pointSize);

            assert(thisModelMatrix == core::Mat4f::eye());

            program->setUniformValue("panoramaCenter", MakeQVec(_projectionCenter));

            program->setUniformValue("projectionMatrix", MakeQMatrix(options.camera.projectionMatrix()));
            program->setUniformValue("modelMatrix", MakeQMatrix(thisModelMatrix));
            program->setUniformValue("viewMatrix", MakeQMatrix(options.camera.viewMatrix()));
            program->setUniformValue("tex", 0);
            program->setUniformValue("bwColor", options.bwColor);
            program->setUniformValue("bwTexColor", options.bwTexColor);

            SetAttributeArrayWithTriMeshVertices(program, "position", _mesh.vertices.front().position);
            SetAttributeArrayWithTriMeshVertices(program, "normal", _mesh.vertices.front().normal);
            SetAttributeArrayWithTriMeshVertices(program, "color", _mesh.vertices.front().color);
            SetAttributeArrayWithTriMeshVertices(program, "texCoord", _mesh.vertices.front().texCoord);
            SetAttributeArrayWithTriMeshVertices(program, "entityIndex", _mesh.vertices.front().entityIndex);
            SetAttributeArrayWithTriMeshVertices(program, "isSelected", _mesh.vertices.front().isSelected);

            program->enableAttributeArray("position");
            program->enableAttributeArray("normal");
            program->enableAttributeArray("color");
            program->enableAttributeArray("texCoord");
            program->enableAttributeArray("entityIndex");
            program->enableAttributeArray("isSelected");

            if (options.renderMode & RenderModeFlag::Triangles) {
                DrawElements(GL_TRIANGLES, _mesh.iTriangles);
            }
            if (options.renderMode & RenderModeFlag::Lines) {
                DrawElements(GL_LINES, _mesh.iLines);
            }
            if (options.renderMode & RenderModeFlag::Points) {
                DrawElements(GL_POINTS, _mesh.iPoints);
            }

            program->disableAttributeArray("position");
            program->disableAttributeArray("normal");
            program->disableAttributeArray("color");
            program->disableAttributeArray("texCoord");
            program->disableAttributeArray("entityIndex");
            program->disableAttributeArray("isSelected");
            
            program->release();

        }

        bool VisualObject::isSingleEntityMeshTriangle(TriMesh::TriangleHandle t) const{
            TriMesh::VertHandle v1, v2, v3;
            _mesh.fetchTriangleVerts(t, v1, v2, v3);
            auto & vs = _mesh.vertices;
            return (vs[v1].entityIndex == vs[v2].entityIndex && vs[v2].entityIndex == vs[v3].entityIndex);
        }

        int VisualObject::entityIDOfMeshTriangle(TriMesh::TriangleHandle t) const {
            TriMesh::VertHandle v1, v2, v3;
            _mesh.fetchTriangleVerts(t, v1, v2, v3);
            auto & vs = _mesh.vertices;
            assert(vs[v1].entityIndex == vs[v2].entityIndex && vs[v2].entityIndex == vs[v3].entityIndex);
            return vs[v1].entityIndex;
        }


        void VisualObject::selectEntity(int entID) {
            _selectedEntities.insert(entID);
            for (auto & v : _mesh.vertices){
                if (v.entityIndex == entID){
                    v.isSelected = true;
                }
            }
        }

        void VisualObject::switchEntitySelection(int entID) {
            if (core::Contains(_selectedEntities, entID)){
                _selectedEntities.erase(entID);
                for (auto & v : _mesh.vertices){
                    if (v.entityIndex == entID){
                        v.isSelected = false;
                    }
                }
            }
            else{
                _selectedEntities.insert(entID);
                for (auto & v : _mesh.vertices){
                    if (v.entityIndex == entID){
                        v.isSelected = true;
                    }
                }
            }
        }

        void VisualObject::clearSelection(){
            _selectedEntities.clear();
            for (auto & v : _mesh.vertices){
                v.isSelected = false;
            }
        }




        class VisualObjectSceneInternal {
            struct VisualObjectMeshTriangleBoundingBoxFunctor {
                inline VisualObjectMeshTriangleBoundingBoxFunctor(VisualObjectSceneInternal * const s) : self(s) {}
                core::Box3 operator()(const VisualObjectMeshTriangle & mti) const {
                    return self->calculatedBoundingBoxesOfMeshTriangles.at(mti);
                }
                VisualObjectSceneInternal * const self;
            };
        public:
            inline VisualObjectSceneInternal()
                : rtree(VisualObjectMeshTriangleBoundingBoxFunctor(this))
            {}

            std::map<VisualObjectHandle, core::Mat4f> calculatedModelMatrices;
            std::map<VisualObjectHandle, core::Box3> calculatedBoundingBoxes;
            std::map<VisualObjectMeshTriangle, core::Box3> calculatedBoundingBoxesOfMeshTriangles;
            core::RTreeWrapper<VisualObjectMeshTriangle, VisualObjectMeshTriangleBoundingBoxFunctor> rtree;
            core::Box3 boundingBox;
        };

        VisualObjectScene::VisualObjectScene() : _internal(new VisualObjectSceneInternal()) {}

        VisualObjectScene::VisualObjectScene(const VisualObjectTree & tree) : _internal(new VisualObjectSceneInternal()) {
            install(tree);
        }

        VisualObjectScene::~VisualObjectScene() {
            delete _internal;
        }

        void VisualObjectScene::update(){            
            clear();
            std::unordered_set<VisualObjectHandle> visited;

            while (true){
                VisualObjectHandle notVisited;
                for (auto & n : _tree.nodes()){
                    if (!Contains(visited, n.topo.hd)){
                        notVisited = n.topo.hd;
                        break;
                    }
                }
                if (notVisited.isInvalid())
                    break;

                _tree.depthFirstSearch(notVisited, [this, &visited](VisualObjectHandle h) -> bool {
                    visited.insert(h);
                    VisualObject * ro = _tree.data(h).get();
                    // update model matrix
                    if (_tree.isRoot(h)) { // is root
                        assert(!core::Contains(_internal->calculatedModelMatrices, h));
                        _internal->calculatedModelMatrices[h] = ro->modelMatrix();
                    }
                    else {
                        assert(core::Contains(_internal->calculatedModelMatrices, _tree.parent(h)));
                        _internal->calculatedModelMatrices[h] = _internal->calculatedModelMatrices.at(_tree.parent(h)) * ro->modelMatrix();
                    }
                    // update bounding box in world space
                    auto & mesh = ro->mesh();
                    std::vector<Point3> transformedVertexPositions;
                    transformedVertexPositions.reserve(mesh.vertices.size());
                    for (auto & vert : mesh.vertices){
                        Point4 transformedHCorner = _internal->calculatedModelMatrices.at(h) * vert.position;
                        transformedVertexPositions.push_back(Point3(transformedHCorner[0], transformedHCorner[1], transformedHCorner[2]) / transformedHCorner[3]);
                    }
                    for (TriMesh::TriangleHandle i = 0; i < mesh.numberOfTriangles(); i++){
                        TriMesh::VertHandle v1, v2, v3;
                        mesh.fetchTriangleVerts(i, v1, v2, v3);
                        _internal->calculatedBoundingBoxesOfMeshTriangles[std::make_pair(h, i)] =
                            BoundingBoxOfContainer(std::vector<Point3>{
                            transformedVertexPositions[v1],
                                transformedVertexPositions[v2],
                                transformedVertexPositions[v3]
                        });
                        _internal->rtree.insert(std::make_pair(h, i));
                    }
                    _internal->calculatedBoundingBoxes[h] = BoundingBoxOfContainer(transformedVertexPositions);
                    return true;
                });

            }

            _internal->boundingBox = 
                core::BoundingBoxOfPairRange(_internal->calculatedBoundingBoxes.begin(),
                _internal->calculatedBoundingBoxes.end());

        }

        void VisualObjectScene::clear() {
            _internal->calculatedModelMatrices.clear();
            _internal->calculatedBoundingBoxesOfMeshTriangles.clear();
            _internal->calculatedBoundingBoxes.clear();
            _internal->rtree.clear();
        }


        const Box3 & VisualObjectScene::boundingBox() const {
            return _internal->boundingBox;
        }

        Box3 VisualObjectScene::boundingBoxOfObject(VisualObjectHandle h) const {
            return _internal->calculatedBoundingBoxes.at(h);
        }

        core::Box3 VisualObjectScene::boundingBoxOfTriangleInObjectMesh(const VisualObjectMeshTriangle & omt) const{
            return _internal->calculatedBoundingBoxesOfMeshTriangles.at(omt);
        }

        void VisualObjectScene::initialize() const {
            for (auto & n : _tree.nodes()){
                if (n.exists){
                    _tree.data(n.topo.hd)->initialize();
                }
            }
        }

        void VisualObjectScene::render(const RenderOptions & options) const {
            for (auto & n : _tree.nodes()){
                if (n.exists){
                    _tree.data(n.topo.hd)->render(options, _internal->calculatedModelMatrices.at(n.topo.hd));
                }
            }
        }

        inline Point3f ToAffinePoint(const Vec4f & p){
            return Point3f(p[0] / p[3], p[1] / p[3], p[2] / p[3]);
        }

        template <class TT>
        bool TriangleIntersection(const Vec<TT, 3> &  V1,  // Triangle vertices
            const Vec<TT, 3> &  V2,
            const Vec<TT, 3> &  V3,
            const Vec<TT, 3> &   O,  //Ray origin
            const Vec<TT, 3> &   D,  //Ray direction

            TT* out,
            TT epsilon) {

            Vec<TT, 3> e1, e2;  //Edge1, Edge2
            Vec<TT, 3> P, Q, T;
            TT det, inv_det, u, v;
            TT t;

            //Find vectors for two edges sharing V1
            e1 = V2 - V1;
            e2 = V3 - V1;
            //Begin calculating determinant - also used to calculate u parameter
            P = D.cross(e2);
            //if determinant is near zero, ray lies in plane of triangle
            det = e1.dot(P);
            //NOT CULLING
            if (det > -epsilon && det < epsilon) return false;
            inv_det = 1.f / det;

            //calculate distance from V1 to ray origin
            T = O - V1;

            //Calculate u parameter and test bound
            u = T.dot(P) * inv_det;
            //The intersection lies outside of the triangle
            if (u < 0.f || u > 1.f) return false;

            //Prepare to test v parameter
            Q = T.cross(e1);

            //Calculate V parameter and test bound
            v = D.dot(Q) * inv_det;
            //The intersection lies outside of the triangle
            if (v < 0.f || u + v  > 1.f) return false;

            t = e2.dot(Q) * inv_det;

            if (t > epsilon) { //ray intersection
                *out = t;
                return true;
            }

            // No hit, no win
            return false;
        }

        static const bool g_UseBruteForce = false;

        VisualObjectMeshTriangle VisualObjectScene::pickOnScreen(const RenderOptions & options,
            const core::Point2 & pOnScreen) const {

            InfiniteLine3 centerRay(options.camera.eye(), normalize(options.camera.spatialDirection(pOnScreen)));

            ///// add ray
            //{
            //    vis::DiscretizeOptions dopts;
            //    dopts.color = vis::ColorFromTag(vis::ColorTag::Blue);
            //    dopts.defaultShaderSource = vis::PredefinedShaderSource(vis::OpenGLShaderSourceDescriptor::XLines);
            //    _tree.add(_tree.firstRoot(), Visualize(Line3(centerRay.anchor, centerRay.anchor + centerRay.direction * 1000), dopts));
            //    update();
            //    initialize();
            //}

            Box3 bboxAll = boundingBox();
            auto bballAll = bboxAll.outerSphere();
            
            // discretize the ray
            double startLen = std::max(0.0, Distance(bballAll.center, options.camera.eye()) - bballAll.radius);
            double stopLen = std::max(0.0, Distance(bballAll.center, options.camera.eye()) + bballAll.radius);
            double stepNum = 1000.0;
            double stepLen = (stopLen - startLen) / stepNum;

            float epsilon = 1e-20f;
            std::map<VisualObjectMeshTriangle, double> resultsWithDepths;

            if (g_UseBruteForce){
                for (auto & mtiBB : _internal->calculatedBoundingBoxesOfMeshTriangles){
                    auto & modelMat = _internal->calculatedModelMatrices.at(mtiBB.first.first);
                    auto & mti = mtiBB.first;
                    auto & vo = _tree.data(mti.first);
                    TriMesh::VertHandle v1, v2, v3;
                    vo->mesh().fetchTriangleVerts(mti.second, v1, v2, v3);
                    float out = 0.0;
                    bool intersected = TriangleIntersection(
                        ToAffinePoint(modelMat * vo->mesh().vertices[v1].position),
                        ToAffinePoint(modelMat * vo->mesh().vertices[v2].position),
                        ToAffinePoint(modelMat * vo->mesh().vertices[v3].position),
                        ConvertTo<float>(centerRay.anchor), ConvertTo<float>(centerRay.direction),
                        &out, epsilon);
                    if (intersected){
                        resultsWithDepths[mti] = out;
                    }
                }
            }
            else{
                for (int i = 0; i < stepNum; i++){
                    double centerLen = startLen + stepLen * i;
                    double nextCenterLen = centerLen + stepLen;

                    Point3 centerP = centerRay.anchor + centerLen * centerRay.direction;
                    Point3 nextCenterP = centerRay.anchor + nextCenterLen * centerRay.direction;

                    Box3 detectionBox = BoundingBox(centerP) | BoundingBox(nextCenterP);
                    detectionBox.expand(stepLen);

                    _internal->rtree.search(detectionBox,
                        [&resultsWithDepths, this, &centerRay, &options, &pOnScreen, epsilon](const VisualObjectMeshTriangle & mti){
                        auto & vo = _tree.data(mti.first);
                        TriMesh::VertHandle v1, v2, v3;
                        vo->mesh().fetchTriangleVerts(mti.second, v1, v2, v3);
                        float out = 0.0;
                        bool intersected = TriangleIntersection(ToAffinePoint(vo->mesh().vertices[v1].position),
                            ToAffinePoint(vo->mesh().vertices[v2].position),
                            ToAffinePoint(vo->mesh().vertices[v3].position),
                            ConvertTo<float>(centerRay.anchor), ConvertTo<float>(centerRay.direction),
                            &out, epsilon);
                        if (intersected){
                            resultsWithDepths[mti] = out;
                        }
                        return true;
                    });
                }
            }

            double minDist = std::numeric_limits<double>::max();
            
            VisualObjectMeshTriangle nearest;

            for (auto & r : resultsWithDepths){
                double d = r.second;
                if (d < minDist){
                    minDist = d;
                    nearest = r.first;
                }
            }

            return nearest;
        }









        // visualizer widget

        class VisualizerWidget : public QGLWidget {
        public:
            RenderOptions options;
            VisualObjectScene scene;

            VisualizerWidget(const Visualizer & v, QWidget * parent = nullptr) 
                : QGLWidget(parent), options(v.renderOptions), scene(v.tree()) {
                setMouseTracking(true);
                setAutoBufferSwap(false);
                grabKeyboard();
            }

        protected:
            void initializeGL() {
                makeCurrent();
                glEnable(GL_MULTISAMPLE);
                GLint bufs;
                GLint samples;
                glGetIntegerv(GL_SAMPLE_BUFFERS, &bufs);
                glGetIntegerv(GL_SAMPLES, &samples);
                qDebug("Have %d buffers and %d samples", bufs, samples);
                qglClearColor(MakeQColor(options.backgroundColor));   
                scene.initialize();
            }

            void paintGL() {
                QPainter painter;
                painter.begin(this);

                painter.beginNativePainting();
                qglClearColor(MakeQColor(options.backgroundColor));

                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                glFrontFace(GL_CW); // face direction set to clockwise
                glEnable(GL_MULTISAMPLE);
                glEnable(GL_DEPTH_TEST);
                glEnable(GL_STENCIL_TEST);

                glEnable(GL_ALPHA_TEST);

                glEnable(GL_BLEND);
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

                glEnable(GL_PROGRAM_POINT_SIZE);
                glEnable(GL_POINT_SPRITE);

                core::PerspectiveCamera & camera = options.camera;
                camera.resizeScreen(core::Size(width(), height()));

                scene.render(options);

                glDisable(GL_DEPTH_TEST);
                glDisable(GL_CULL_FACE);

                painter.endNativePainting();
                swapBuffers();
            }

            void resizeGL(int w, int h) {
                core::PerspectiveCamera & camera = options.camera;
                camera.resizeScreen(core::Size(w, h));
                glViewport(0, 0, w, h);
            }

        public:
            void autoSetCamera() {
                auto sphere = scene.boundingBox().outerSphere();
                options.camera.resizeScreen(core::Size(width(), height()), false);
                options.camera.focusOn(sphere, true);
                update();
            }

        protected:
            virtual void mousePressEvent(QMouseEvent * e) override {
                _lastPos = e->pos();
                if (e->buttons() & Qt::RightButton)
                    setCursor(Qt::OpenHandCursor);
                else if (e->buttons() & Qt::MidButton)
                    setCursor(Qt::SizeAllCursor);
                else if (e->buttons() & Qt::LeftButton){
                    VisualObjectHandle oh;
                    TriMesh::TriangleHandle t;
                    std::tie(oh, t) = scene.pickOnScreen(options, core::Point2(e->pos().x(), e->pos().y()));                    
                    if (oh.isValid()){
                        int entityID = scene.tree().data(oh)->entityIDOfMeshTriangle(t);
                        if (e->modifiers() & Qt::ControlModifier){
                            scene.switchSelect(std::make_pair(oh, entityID));
                        }
                        else{
                            scene.clearSelection();
                            scene.select(std::make_pair(oh, entityID));
                        }
                        scene.tree().data(oh)
                            ->invokeCallbackFunction(InteractionID::ClickLeftButton, scene.tree(), std::make_pair(oh, entityID));
                    }
                    else if (!(e->modifiers() & Qt::ControlModifier)){
                        scene.clearSelection();
                    }
                }
                update();
            }

            virtual void mouseMoveEvent(QMouseEvent * e) override {
                QVector3D t(e->pos() - _lastPos);
                t.setX(-t.x());
                auto sphere = scene.boundingBox().outerSphere();
                if (e->buttons() & Qt::RightButton) {
                    core::Vec3 trans = t.x() * options.camera.rightward() + t.y() * options.camera.upward();
                    trans *= 0.02;
                    options.camera.moveEyeWithCenterFixed(trans, sphere, true, true);
                    setCursor(Qt::ClosedHandCursor);
                    update();
                }
                else if (e->buttons() & Qt::MidButton) {
                    core::Vec3 trans = t.x() * options.camera.rightward() + t.y() * options.camera.upward();
                    trans *= 0.02;
                    options.camera.translate(trans, sphere, true);
                    update();
                }
                _lastPos = e->pos();
            }

            virtual void wheelEvent(QWheelEvent * e) override {
                auto sphere = scene.boundingBox().outerSphere();
                double d = e->delta() / 10;
                double dist = core::Distance(options.camera.eye(), options.camera.center());
                core::Vec3 trans = d * dist / 1000.0 * options.camera.forward();
                options.camera.moveEyeWithCenterFixed(trans, sphere, false, true);
                update();
            }

            virtual void mouseReleaseEvent(QMouseEvent * e) override {
                unsetCursor();
            }

            virtual void keyPressEvent(QKeyEvent * e) override {
                if (e->key() == Qt::Key_Space){
                    for (auto & n : scene.tree().nodes()){
                        if (n.exists){
                            for (int entityID : n.data->selectedEntities()){
                                n.data->invokeCallbackFunction(InteractionID::PressSpace,
                                    scene.tree(),
                                    VisualObjectEntityID{ n.topo.hd, entityID });
                            }
                        }
                    }
                }
            }

        private:
            QPointF _lastPos;
        };


        class VisualizerMainWindow : public QMainWindow {
        public:
            explicit VisualizerMainWindow(QWidget * parent = nullptr) : QMainWindow(parent) {
                setupGui();
            }
            void setupGui() {
                auto menuView = this->menuBar()->addMenu(tr("View"));
                auto menuAbout = this->menuBar()->addMenu(tr("About"));
                this->statusBar()->show();

                auto actionAbout = menuAbout->addAction(tr("About"));
                connect(actionAbout, &QAction::triggered, [this](){
                    QMessageBox::about(this, tr("About this program"),
                        tr("Panoramix.Vis is the visulization module of project Panoramix developped by Yang Hao."));
                });
            }
        };

        







        void Visualizer::show(bool doModal, bool autoSetCamera) {
            auto app = Singleton::InitGui();
            VisualizerWidget * w = new VisualizerWidget(*this);
            VisualizerMainWindow * mwin = new VisualizerMainWindow();
            mwin->setCentralWidget(w);
            mwin->setAttribute(Qt::WA_DeleteOnClose);
            mwin->resize(MakeQSize(renderOptions.camera.screenSize()));
            mwin->setWindowTitle(QString::fromStdString(renderOptions.winName));
            mwin->setWindowIcon(Singleton::DefaultConfiguration().icon);
            mwin->setStyleSheet(Singleton::DefaultConfiguration().css);
            auto palette = mwin->palette();
            palette.setColor(QPalette::Window, MakeQColor(renderOptions.backgroundColor));
            mwin->setPalette(palette);
            //qDebug() << mwin->styleSheet();
            if (autoSetCamera) {
                w->autoSetCamera();
            }
            mwin->show();
            if (doModal) {
                Singleton::ContinueGui(); // qApp->exec()
            }
        }

    }
}