#include <QtOpenGL>

#include "../core/containers.hpp"
#include "qt_glue.hpp"
#include "scene.hpp"

namespace panoramix {
    namespace gui {

        using namespace core;


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
                program(nullptr) {}
            inline void initialize() { if(!program) program = new QOpenGLShaderProgram; }
            inline bool isLinked() const { return program ? program->isLinked() : false; }
            inline QString log() const { return program ? program->log() : "QOpenGLShaderProgram not initialized!"; }
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
                    qDebug() << (_internal->log());
            }

            delete _internal;
        }

        void VisualObject::setShaderSource(const OpenGLShaderSource & shaderSource){
            if (_internal->isLinked()){
                qDebug() << "program is already linked! setting shaders failed!";
                return;
            }
            _shaderSource = shaderSource;
        }

        void VisualObject::initialize() const {

            auto vo = _internal;
            vo->initialize();
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
            vo->initialize();
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




        class SceneInternal {
            struct VisualObjectMeshTriangleBoundingBoxFunctor {
                inline VisualObjectMeshTriangleBoundingBoxFunctor(SceneInternal * const s) : self(s) {}
                core::Box3 operator()(const VisualObjectMeshTriangle & mti) const {
                    return self->calculatedBoundingBoxesOfMeshTriangles.at(mti);
                }
                SceneInternal * const self;
            };
        public:
            inline SceneInternal()
                : rtree(VisualObjectMeshTriangleBoundingBoxFunctor(this))
            {}

            std::map<VisualObjectHandle, core::Mat4f> calculatedModelMatrices;
            std::map<VisualObjectHandle, core::Box3> calculatedBoundingBoxes;
            std::map<VisualObjectMeshTriangle, core::Box3> calculatedBoundingBoxesOfMeshTriangles;
            core::RTreeWrapper<VisualObjectMeshTriangle, VisualObjectMeshTriangleBoundingBoxFunctor> rtree;
            core::Box3 boundingBox;
        };

        Scene::Scene() : _internal(new SceneInternal()) {}

        Scene::Scene(const VisualObjectTree & tree) : _internal(new SceneInternal()) {
            install(tree);
        }

        Scene::~Scene() {
            delete _internal;
        }

        void Scene::update(){
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
                if (notVisited.invalid())
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

        void Scene::clear() {
            _internal->calculatedModelMatrices.clear();
            _internal->calculatedBoundingBoxesOfMeshTriangles.clear();
            _internal->calculatedBoundingBoxes.clear();
            _internal->rtree.clear();
        }


        const Box3 & Scene::boundingBox() const {
            return _internal->boundingBox;
        }

        Box3 Scene::boundingBoxOfObject(VisualObjectHandle h) const {
            return _internal->calculatedBoundingBoxes.at(h);
        }

        core::Box3 Scene::boundingBoxOfTriangleInObjectMesh(const VisualObjectMeshTriangle & omt) const{
            return _internal->calculatedBoundingBoxesOfMeshTriangles.at(omt);
        }

        void Scene::initialize() const {
            for (auto & n : _tree.nodes()){
                if (n.exists){
                    _tree.data(n.topo.hd)->initialize();
                }
            }
        }

        void Scene::render(const RenderOptions & options) const {
            glFrontFace(GL_CCW); // face direction set to clockwise
            glEnable(GL_MULTISAMPLE);
            glEnable(GL_DEPTH_TEST);
            glEnable(GL_STENCIL_TEST);

            glEnable(GL_ALPHA_TEST);

            if (options.cullFrontFace || options.cullBackFace){
                glEnable(GL_CULL_FACE);
                if (options.cullFrontFace && !options.cullBackFace)
                    glCullFace(GL_FRONT);
                else if (!options.cullFrontFace && options.cullBackFace)
                    glCullFace(GL_BACK);
                else
                    glCullFace(GL_FRONT_AND_BACK);
            }
            else{
                glDisable(GL_CULL_FACE);
            }

            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

            glEnable(GL_PROGRAM_POINT_SIZE);
            glEnable(GL_POINT_SPRITE);

            for (auto & n : _tree.nodes()){
                if (n.exists){
                    _tree.data(n.topo.hd)->render(options, _internal->calculatedModelMatrices.at(n.topo.hd));
                }
            }

            glDisable(GL_DEPTH_TEST);
            glDisable(GL_CULL_FACE);
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

        VisualObjectMeshTriangle Scene::pickOnScreen(const RenderOptions & options,
            const core::Point2 & pOnScreen) const {

            Ray3 centerRay(options.camera.eye(), normalize(options.camera.spatialDirection(pOnScreen)));

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
                        vec_cast<float>(centerRay.anchor), vec_cast<float>(centerRay.direction),
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
                            vec_cast<float>(centerRay.anchor), vec_cast<float>(centerRay.direction),
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




        PerspectiveCamera Scene::perfectView(int width, int height, const Vec3 & up) const{
            PerspectiveCamera camera;
            auto sphere = boundingBox().outerSphere();
            camera.setUp(up, false);
            camera.resizeScreen(core::Size(width, height), false);
            camera.focusOn(sphere, true);
            return camera;
        }

 
    }
}