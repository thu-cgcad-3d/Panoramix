#include "qt_opengl_object.hpp"

namespace panoramix {
    namespace vis {

        // opengl _mesh data implementation
        OpenGLMeshData::VertHandle OpenGLMeshData::addVertex(const OpenGLMeshData::Vertex & v) {
            vertices.push_back(v);
            iPoints.push_back(vertices.size() - 1);
            return vertices.size() - 1;
        }

        OpenGLMeshData::VertHandle OpenGLMeshData::addVertex(const QVector4D & p, const QVector3D & n, const QVector4D & c, const QVector2D & t) {
            Vertex v;
            v.position4 = p;
            v.normal3 = n;
            v.color4 = c;
            v.texCoord2 = t;
            return addVertex(v);
        }

        OpenGLMeshData::LineHandle OpenGLMeshData::addLine(OpenGLMeshData::VertHandle v1, OpenGLMeshData::VertHandle v2) {
            iLines.push_back(v1);
            iLines.push_back(v2);
            return iLines.size() / 2;
        }

        OpenGLMeshData::LineHandle OpenGLMeshData::addIsolatedLine(const Vertex & v1, const Vertex & v2) {
            vertices.push_back(v1);
            iLines.push_back(vertices.size() - 1);
            vertices.push_back(v2);
            iLines.push_back(vertices.size() - 1);
            return iLines.size() / 2;
        }

        OpenGLMeshData::TriangleHandle OpenGLMeshData::addTriangle(OpenGLMeshData::VertHandle v1, OpenGLMeshData::VertHandle v2, OpenGLMeshData::VertHandle v3) {
            iTriangles.push_back(v1);
            iTriangles.push_back(v2);
            iTriangles.push_back(v3);
            return iTriangles.size() / 3;
        }

        OpenGLMeshData::TriangleHandle OpenGLMeshData::addIsolatedTriangle(const Vertex & v1, const Vertex & v2, const Vertex & v3) {
            vertices.push_back(v1);
            iTriangles.push_back(vertices.size() - 1);
            vertices.push_back(v2);
            iTriangles.push_back(vertices.size() - 1);
            vertices.push_back(v3);
            iTriangles.push_back(vertices.size() - 1);
            return iTriangles.size() / 3;
        }

        void OpenGLMeshData::addQuad(OpenGLMeshData::VertHandle v1, OpenGLMeshData::VertHandle v2, OpenGLMeshData::VertHandle v3, OpenGLMeshData::VertHandle v4){
            addTriangle(v1, v2, v3);
            addTriangle(v1, v3, v4);
        }

        namespace {
            // algorithms
            inline double tDet(double* data) {
                double tmp1 = data[0 * 3 + 0] * (data[1 * 3 + 1] * data[2 * 3 + 2] - data[1 * 3 + 2] * data[2 * 3 + 1]);
                double tmp2 = data[0 * 3 + 1] * (data[1 * 3 + 0] * data[2 * 3 + 2] - data[1 * 3 + 2] * data[2 * 3 + 0]);
                double tmp3 = data[0 * 3 + 2] * (data[1 * 3 + 0] * data[2 * 3 + 1] - data[1 * 3 + 1] * data[2 * 3 + 0]);
                return tmp1 - tmp2 + tmp3;
            }

            template <typename PointT>
            inline bool tLeft(const PointT& p, const PointT& a, const PointT& b) {
                double data[9] = { a[0], a[1], 1, b[0], b[1], 1, p[0], p[1], 1 };
                return tDet(data) > 0;
            }

            template <typename PointT>
            inline bool tInTriangle(const PointT& p, const PointT& a, const PointT& b, const PointT& c) {
                bool lab = tLeft(p, a, b);
                bool lbc = tLeft(p, b, c);
                bool lca = tLeft(p, c, a);
                return lab == lbc && lbc == lca;
            }

            template <typename PointT>
            inline double tSqDist(const PointT& p1, const PointT& p2) {
                auto sub = p1 - p2;
                return sub[0] * sub[0] + sub[1] * sub[1] + sub[2] * sub[2];
            }

#define qRoundNear(a, b, size) \
    (abs(a - b) <= 1 || a == 0 && b == (size)-1 || a == (size)-1 && b == 0)

            template <typename VHandleT, typename VHandleGetPointFunctorT>
            QList<VHandleT> tTriangulate(VHandleGetPointFunctorT _mesh, const QList<VHandleT>& vhs) {
                QList<VHandleT> triangles;
                QQueue<QList<int> > vhIndexGroupQ;
                QList<int> indexG;
                for (int i = 0; i < vhs.size(); i++)
                    indexG.push_back(i);
                vhIndexGroupQ.push_back(indexG);

                while (!vhIndexGroupQ.empty()) {
                    QList<int> is = vhIndexGroupQ.first();
                    vhIndexGroupQ.pop_front();

                    Q_ASSERT(is.size() >= 3);
                    if (is.size() <= 2)
                        continue;

                    if (is.size() == 3)
                        triangles << vhs[is[0]] << vhs[is[1]] << vhs[is[2]];
                    else{
                        // leftmost
                        int leftmostII = 0;
                        auto leftmostP = _mesh(vhs[is[leftmostII]]);
                        for (int i = 0; i < is.size(); i++){
                            auto p = _mesh(vhs[is[i]]);
                            if (p[0] < leftmostP[0]){
                                leftmostII = i;
                                leftmostP = p;
                            }
                        }

                        int leftmostPrevII = (leftmostII + is.size() - 1) % is.size();
                        int leftmostNextII = (leftmostII + 1) % is.size();
                        auto a = _mesh(vhs[is[leftmostPrevII]]);
                        auto b = _mesh(vhs[is[leftmostNextII]]);

                        int innerLeftmostII = -1;
                        decltype(a) innerLeftmostP;
                        for (int i = 0; i < is.size(); i++){
                            if (qRoundNear(i, leftmostII, is.size()))
                                continue;
                            auto p = _mesh(vhs[is[i]]);
                            if (tInTriangle(p, a, leftmostP, b))
                            {
                                if (innerLeftmostII == -1){
                                    innerLeftmostII = i;
                                    innerLeftmostP = p;
                                }
                                else if (p[0] < innerLeftmostP[0]){
                                    innerLeftmostII = i;
                                    innerLeftmostP = p;
                                }
                            }
                        }

                        int split1 = leftmostII;
                        int split2 = innerLeftmostII;
                        if (innerLeftmostII < 0){
                            split1 = leftmostPrevII;
                            split2 = leftmostNextII;
                        }

                        Q_ASSERT(split1 != split2);

                        QList<int> part1, part2;

                        for (int i = split1; i != split2; i = (i + 1) % is.size())
                            part1.push_back(is[i]);
                        part1.push_back(is[split2]);
                        for (int i = split2; i != split1; i = (i + 1) % is.size())
                            part2.push_back(is[i]);
                        part2.push_back(is[split1]);

                        Q_ASSERT(part1.size() >= 3);
                        Q_ASSERT(part2.size() >= 3);

                        is.clear();

                        vhIndexGroupQ.push_back(part1);
                        vhIndexGroupQ.push_back(part2);
                    }
                }

                return triangles;
            }
        }

        void OpenGLMeshData::addPolygon(const QList<OpenGLMeshData::VertHandle> & vhs) {
            Q_ASSERT(vhs.size() >= 3);
            // get normal direction
            QVector3D normal = QVector3D::crossProduct(
                (vertices[vhs[1]].position4.toVector3DAffine() - vertices[vhs[0]].position4.toVector3DAffine()),
                (vertices[vhs[2]].position4.toVector3DAffine() - vertices[vhs[1]].position4.toVector3DAffine()));
            normal.normalize();

            auto pointGetter = [&](VertHandle vh){
                QVector3D v = vertices[vh].position4.toVector3DAffine();
                return (v - QVector3D::dotProduct(v, normal) * normal).toVector2D();
            };

            QList<VertHandle> triangles = tTriangulate(pointGetter, vhs);
            for (int i = 0; i < triangles.size(); i += 3){
                addTriangle(triangles[i], triangles[i + 1], triangles[i + 2]);
            }
        }

        void OpenGLMeshData::clear(){
            vertices.clear();
            iPoints.clear();
            iLines.clear();
            iTriangles.clear();
        }

        QPair<QVector3D, QVector3D> OpenGLMeshData::boundingBox() const {
            if (vertices.empty())
                return QPair<QVector3D, QVector3D>();
            QPair<QVector3D, QVector3D> box(vertices.front().position4.toVector3DAffine(), 
                vertices.front().position4.toVector3DAffine());
            for (auto & v : vertices){
                auto p = v.position4.toVector3DAffine();
                for (int i = 0; i < 3; i++){
                    if (box.first[i] > p[i])
                        box.first[i] = p[i];
                    if (box.second[i] < p[i])
                        box.second[i] = p[i];
                }
            }
            return box;
        }


        namespace {

            static const OpenGLShaderSource SSNormalPoints = {
                "attribute highp vec4 position;\n"
                "attribute lowp float pointSize;\n"
                "attribute highp vec3 normal;\n"
                "attribute lowp vec4 color;\n"
                "attribute lowp vec2 texCoord;\n"
                "uniform highp mat4 viewMatrix;\n"
                "uniform highp mat4 modelMatrix;\n"
                "uniform highp mat4 projectionMatrix;\n"
                "varying vec4 pixelColor;\n"
                "void main(void)\n"
                "{\n"
                "    gl_Position = projectionMatrix * viewMatrix * modelMatrix * position;\n"
                "    gl_PointSize = pointSize;\n"
                "    pixelColor = color;\n"
                "}\n",

                "varying lowp vec4 pixelColor;\n"
                "void main(void)\n"
                "{\n"
                "   gl_FragColor = pixelColor;\n"
                "   lowp float distance = length(gl_PointCoord - vec2(0.5));\n"
                "   if(distance > 0.4 && distance <= 0.5)\n"
                "       gl_FragColor.a = 1.0 - (distance - 0.4) * 0.1;\n"
                "   else if(distance > 0.5)\n"
                "       discard;\n"
                "}\n"

                ""
            };

            static const OpenGLShaderSource SSNormalLines = {
                "attribute highp vec4 position;\n"
                "attribute lowp float pointSize;\n"
                "attribute highp vec3 normal;\n"
                "attribute lowp vec4 color;\n"
                "attribute lowp vec2 texCoord;\n"
                "uniform highp mat4 viewMatrix;\n"
                "uniform highp mat4 modelMatrix;\n"
                "uniform highp mat4 projectionMatrix;\n"
                "varying vec4 pixelColor;\n"
                "void main(void)\n"
                "{\n"
                "    gl_Position = projectionMatrix * viewMatrix * modelMatrix * position;\n"
                "    pixelColor = color;\n"
                "}\n",

                "varying lowp vec4 pixelColor;\n"
                "void main(void)\n"
                "{\n"
                "    gl_FragColor = pixelColor;\n"
                "}\n"
            };

            static const OpenGLShaderSource SSNormalTriangles = {
                "attribute highp vec4 position;\n"
                "attribute lowp float pointSize;\n"
                "attribute highp vec3 normal;\n"
                "attribute lowp vec4 color;\n"
                "attribute lowp vec2 texCoord;\n"
                "uniform highp mat4 viewMatrix;\n"
                "uniform highp mat4 modelMatrix;\n"
                "uniform highp mat4 projectionMatrix;\n"
                "varying vec4 pixelColor;\n"
                "void main(void)\n"
                "{\n"
                "    gl_Position = projectionMatrix * viewMatrix * modelMatrix * position;\n"
                "    pixelColor = color;\n"
                "}\n",

                "varying lowp vec4 pixelColor;\n"
                "void main(void)\n"
                "{\n"
                "    gl_FragColor = pixelColor;\n"
                "}\n"
            };

        }



        OpenGLShaderSource PredefinedShaderSource(OpenGLShaderSourceName name) {
            switch (name) {
            case OpenGLShaderSourceName::NormalPoints: return SSNormalPoints;
            case OpenGLShaderSourceName::NormalLines: return SSNormalLines;
            case OpenGLShaderSourceName::NormalTriangles: return SSNormalTriangles;
            default:
                return OpenGLShaderSource();
            }
        }





        // opengl object implementation
        OpenGLObject::OpenGLObject(QObject *parent) : QObject(parent), _texture(NULL) {
            _program = new QOpenGLShaderProgram(this);
            _vertexArrayBuffer = _pointsIndicesBuffer = _linesIndicesBuffer = _trianglesIndicesBuffer = 0;
            glGenBuffers(1, &_vertexArrayBuffer);
            glGenBuffers(1, &_pointsIndicesBuffer);
            glGenBuffers(1, &_linesIndicesBuffer);
            glGenBuffers(1, &_trianglesIndicesBuffer);
        }

        OpenGLObject::~OpenGLObject() {
            delete _texture;
            GLuint buffers[] = { _vertexArrayBuffer, _pointsIndicesBuffer, _linesIndicesBuffer, _trianglesIndicesBuffer };
            glDeleteBuffers(4, buffers);
        }

        void OpenGLObject::setUpShaders(const OpenGLShaderSource & ss){
            if (!_program->addShaderFromSourceCode(QOpenGLShader::Vertex, ss.vertexShaderSource))
                error(_program->log());
            if (!_program->addShaderFromSourceCode(QOpenGLShader::Fragment, ss.fragmentShaderSource))
                error(_program->log());
            /*if (!_program->addShaderFromSourceCode(QOpenGLShader::Geometry, ss.geometryShaderSource))
                error(_program->log());*/
            if (!_program->link() || !_program->bind()){
                error(_program->log());
                return;
            }

            Q_ASSERT(_program->isLinked());
            qDebug() << _program->log();
            _program->release();
        }


        void OpenGLObject::setUpTexture(const QImage & tex) {
            _program->bind();
            Q_ASSERT(_program->isLinked());
            if (_texture){
                delete _texture;
            }
            _texture = new QOpenGLTexture(tex.mirrored());
            _texture->bind();
            _texture->setMinificationFilter(QOpenGLTexture::Linear);
            _texture->setMagnificationFilter(QOpenGLTexture::Linear);
            _texture->release();
            _program->release();
        }

        void OpenGLObject::setUpMesh(const OpenGLMeshData & m){
            _mesh = m;

            glBindBuffer(GL_ARRAY_BUFFER, _vertexArrayBuffer);
            glBufferData(GL_ARRAY_BUFFER, _mesh.vertices.size() * sizeof(OpenGLMeshData::Vertex), _mesh.vertices.data(), GL_STATIC_DRAW);            
           
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _pointsIndicesBuffer);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, _mesh.iPoints.size() * sizeof(OpenGLMeshData::VertHandle), _mesh.iPoints.data(), GL_STATIC_DRAW);
            
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _linesIndicesBuffer);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, _mesh.iLines.size() * sizeof(OpenGLMeshData::VertHandle), _mesh.iLines.data(), GL_STATIC_DRAW);
            
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _trianglesIndicesBuffer);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, _mesh.iTriangles.size() * sizeof(OpenGLMeshData::VertHandle), _mesh.iTriangles.data(), GL_STATIC_DRAW);

        }

        void OpenGLObject::render(RenderModeFlags mode, const QMatrix4x4 & projection,
            const QMatrix4x4 & view, const QMatrix4x4 & model) {

            

            if (_mesh.vertices.isEmpty())
                return;

            Q_ASSERT(_program->isLinked());
            _program->bind();

            if (_texture && _texture->isCreated())
                _texture->bind(0);

            glLineWidth(_mesh.vertices.first().lineWidth1);
            glBindBuffer(GL_ARRAY_BUFFER, _vertexArrayBuffer);

            _program->setUniformValue("projectionMatrix", projection);
            _program->setUniformValue("viewMatrix", view);
            _program->setUniformValue("modelMatrix", model);
            _program->setUniformValue("tex", 0);
            _program->setUniformValue("panoramaCenter", QVector3D(0, 0, 0));

            glVertexAttribPointer(_program->attributeLocation("position"), 4, GL_FLOAT, GL_FALSE,
                sizeof(OpenGLMeshData::Vertex),
                //OFFSET_OF(_mesh.vertices.first().position4, _mesh.vertices.first().position4));
                (const void*)offsetof(OpenGLMeshData::Vertex, position4));
            glVertexAttribPointer(_program->attributeLocation("normal"), 3, GL_FLOAT, GL_FALSE,
                sizeof(OpenGLMeshData::Vertex),
                //OFFSET_OF(_mesh.vertices.first().normal3, _mesh.vertices.first().position4));
                (const void*)offsetof(OpenGLMeshData::Vertex, normal3));
            glVertexAttribPointer(_program->attributeLocation("color"), 4, GL_FLOAT, GL_FALSE,
                sizeof(OpenGLMeshData::Vertex),
                //OFFSET_OF(_mesh.vertices.first().color4, _mesh.vertices.first().position4));
                (const void*)offsetof(OpenGLMeshData::Vertex, color4));
            glVertexAttribPointer(_program->attributeLocation("texCoord"), 2, GL_FLOAT, GL_FALSE,
                sizeof(OpenGLMeshData::Vertex),
                //OFFSET_OF(_mesh.vertices.first().texCoord2, _mesh.vertices.first().position4));
                (const void*)offsetof(OpenGLMeshData::Vertex, texCoord2));
            glVertexAttribPointer(_program->attributeLocation("pointSize"), 1, GL_FLOAT, GL_FALSE,
                sizeof(OpenGLMeshData::Vertex),
                //OFFSET_OF(_mesh.vertices.first().pointSize1, _mesh.vertices.first().position4));
                (const void*)offsetof(OpenGLMeshData::Vertex, pointSize1));
            glVertexAttribPointer(_program->attributeLocation("lineWidth"), 1, GL_FLOAT, GL_FALSE,
                sizeof(OpenGLMeshData::Vertex),
                //OFFSET_OF(_mesh.vertices.first().lineWidth1, _mesh.vertices.first().position4));
                (const void*)offsetof(OpenGLMeshData::Vertex, lineWidth1));

            _program->enableAttributeArray("position");
            _program->enableAttributeArray("normal");
            _program->enableAttributeArray("color");
            _program->enableAttributeArray("texCoord");
            _program->enableAttributeArray("pointSize");
            _program->enableAttributeArray("lineWidth");

            if (mode & RenderModeFlag::Triangles){
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _trianglesIndicesBuffer);
                glDrawElements(GL_TRIANGLES, _mesh.iTriangles.size(), GL_UNSIGNED_INT, 0);
            }

            if (mode & RenderModeFlag::Lines) {
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _linesIndicesBuffer);
                glDrawElements(GL_LINES, _mesh.iLines.size(), GL_UNSIGNED_INT, 0);
            }

            if (mode & RenderModeFlag::Points){
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _pointsIndicesBuffer);
                glDrawElements(GL_POINTS, _mesh.iPoints.size(), GL_UNSIGNED_INT, 0);
            }

            _program->disableAttributeArray("position");
            _program->disableAttributeArray("normal");
            _program->disableAttributeArray("color");
            _program->disableAttributeArray("texCoord");
            _program->disableAttributeArray("pointSize");
            _program->disableAttributeArray("lineWidth");

            _program->release();
        }

        void OpenGLObject::error(const QString & message) {
            qWarning() << message;
            emit errorOccored(message);
        }


    }
}