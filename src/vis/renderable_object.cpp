#include <QtOpenGL>

#include "../core/misc.hpp"
#include "../core/utilities.hpp"
#include "qt_glue.hpp"

#include "renderable_object.hpp"

#define GL_ALPHA_TEST 0x0BC0

namespace panoramix {
    namespace vis {

        using namespace core;

        RenderableObject::RenderableObject(RenderableObject * parent) : _parent(parent) {
            if (_parent)
                _parent->_children.push_back(this);
            _modelMat = Mat4::eye();
        }

        RenderableObject::~RenderableObject() {
            for (RenderableObject * ch : _children)
                delete ch;
        }

        void RenderableObject::renderWithCamera(RenderModeFlags mode, const core::PerspectiveCamera & cam) const {
            render(mode, cam.viewProjectionMatrix() * _modelMat);
        }




        DefaultRenderState::DefaultRenderState() : 
            foregroundColor(0, 0, 0),
            pointSize(10.0f),
            lineWidth(2.0f),
            colorTable(ColorTableDescriptor::AllColors) 
        {}


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

            template <class T, int N>
            inline void SetAttributeArrayWithOpenGLMeshVertices(
                QOpenGLShaderProgram * program, const char * attributeName,
                const core::Vec<T, N> & firstVector) {
                program->setAttributeArray(attributeName, OpenGLDataTraits<T>::GLType, 
                    &firstVector, N, sizeof(OpenGLMesh::Vertex));
            }

            template <class T>
            inline void DrawElements(GLenum mode, const std::vector<T> & indices) {
                glDrawElements(mode, indices.size(), OpenGLDataTraits<T>::GLType, indices.data());
            }

            class GLObject : public RenderableObject {
            public:
                explicit GLObject(const OpenGLMesh & m, float pointSize, float lineWidth,
                    const OpenGLShaderSource & ss, const core::Vec3 & panoramaCenter, 
                    RenderableObject * parent)
                    : RenderableObject(parent), _mesh(m), _pointSize(pointSize), _lineWidth(lineWidth),
                    _shaderSource(ss), _panoramaCenter(panoramaCenter) {
                    _program = new QOpenGLShaderProgram;
                    _texture = new QOpenGLTexture(QOpenGLTexture::Target2D);
                }

                virtual ~GLObject() {
                    delete _program;
                    delete _texture;
                }

            public:
                virtual void initialize() const override {
                    // setup shaders
                    if (!_program->addShaderFromSourceCode(QOpenGLShader::Vertex,
                        QByteArray::fromRawData(_shaderSource.vertexShaderSource().data(), _shaderSource.vertexShaderSource().size())))
                        qDebug() << (_program->log());
                    if (!_program->addShaderFromSourceCode(QOpenGLShader::Fragment,
                        QByteArray::fromRawData(_shaderSource.fragmentShaderSource().data(), _shaderSource.fragmentShaderSource().size())))
                        qDebug() << (_program->log());
                    if (!_program->link() || !_program->bind()) {
                        qDebug() << (_program->log());
                        return;
                    }

                    Q_ASSERT(_program->isLinked());

                    // setup texture
                    if (!_textureImage.isNull()){
                        if (!_texture->isCreated()){
                            _texture->create();
                        }
                        if (!_texture->isCreated()){
                            qDebug() << _program->log();
                        }
                        Q_ASSERT(_texture->textureId());
                        _texture->setData(_textureImage.mirrored());
                        _texture->setMinificationFilter(QOpenGLTexture::Linear);
                        _texture->setMagnificationFilter(QOpenGLTexture::Linear);
                        _texture->release();
                    }

                    _program->release();
                }

                virtual void render(RenderModeFlags mode, const Mat4 & mat) const override {
                    if (_mesh.vertices().empty())
                        return;

                    Q_ASSERT(_program->isLinked());
                    _program->bind();

                    if (_texture && _texture->isCreated())
                        _texture->bind(0);

                    glLineWidth(_lineWidth);

                    _program->setUniformValue("matrix", MakeQMatrix(mat));
                    _program->setUniformValue("tex", 0);
                    _program->setUniformValue("panoramaCenter", MakeQVec(_panoramaCenter));
                    _program->setUniformValue("pointSize", _pointSize);

                    SetAttributeArrayWithOpenGLMeshVertices(_program, "position", _mesh.vertices().front().position4);
                    SetAttributeArrayWithOpenGLMeshVertices(_program, "normal", _mesh.vertices().front().normal3);
                    SetAttributeArrayWithOpenGLMeshVertices(_program, "color", _mesh.vertices().front().color4);
                    SetAttributeArrayWithOpenGLMeshVertices(_program, "texCoord", _mesh.vertices().front().texCoord2);

                    _program->enableAttributeArray("position");
                    _program->enableAttributeArray("normal");
                    _program->enableAttributeArray("color");
                    _program->enableAttributeArray("texCoord");

                    if (mode & RenderModeFlag::Triangles) {
                        DrawElements(GL_TRIANGLES, _mesh.iTriangles());
                    }
                    if (mode & RenderModeFlag::Lines) {
                        DrawElements(GL_LINES, _mesh.iLines());
                    }
                    if (mode & RenderModeFlag::Points) {
                        DrawElements(GL_POINTS, _mesh.iPoints());
                    }

                    _program->disableAttributeArray("position");
                    _program->disableAttributeArray("normal");
                    _program->disableAttributeArray("color");
                    _program->disableAttributeArray("texCoord");

                    _program->release();
                    
                }

                virtual core::Box3 primaryBoundingBox() const { return BoundingBox(_mesh); }
                virtual float distanceTo(const core::InfiniteLine3 & ray) const {
                    // TODO
                    NOT_IMPLEMENTED_YET();
                }

                virtual void setTexture(const core::Image & im) {
                    _textureImage = vis::MakeQImage(im);
                }

            protected:
                OpenGLMesh _mesh;
                QImage _textureImage;
                OpenGLShaderSource _shaderSource;
                QOpenGLShaderProgram * _program;
                QOpenGLTexture * _texture;
                float _lineWidth;
                float _pointSize;
                core::Vec3 _panoramaCenter;
            };
        }




        namespace {

            class GLPointsObject : public GLObject {
            public:
                template <class PointsIteratorT>
                explicit GLPointsObject(PointsIteratorT && begin, PointsIteratorT && end,
                    float pointSize, const Color & color, RenderableObject * parent)
                    : GLObject(OpenGLMesh::FromPoints(begin, end), pointSize, 1.0f, 
                    OpenGLShaderSourceDescriptor::DefaultPoints, core::Vec3(0, 0, 0), parent), 
                    _points(begin, end), _color(color) {
                    // update color
                    for (auto & v : _mesh.vertices()) {
                        QColor c = MakeQColor(color);
                        v.color4 = core::Vec4f(c.redF(), c.greenF(), c.blueF(), 1.0);
                    }
                }

                virtual float distanceTo(const core::InfiniteLine3 & ray) const {
                    float distance = std::numeric_limits<float>::max();
                    for (auto & p : _points) {
                        float d = core::DistanceFromPointToLine(p, ray).first;
                        if (d < distance)
                            distance = d;
                    }
                    return distance;
                }

            private:
                std::vector<Point3> _points;
                Color _color;
            };
        }

        // point
        RenderableObject * MakeRenderable(const Point3 & p,
            const DefaultRenderState & param, RenderableObject * parent) {
            return new GLPointsObject(&p, &p + 1, param.pointSize, param.foregroundColor, parent);
        }

        // points
        RenderableObject * MakeRenderable(const std::vector<Point3> & points,
            const DefaultRenderState & param, RenderableObject * parent) {
            return new GLPointsObject(points.begin(), points.end(), param.pointSize, param.foregroundColor, parent);
        }




        namespace {

            class GLLinesObject : public GLObject {
            public:
                template <class LinesIteratorT>
                explicit GLLinesObject(LinesIteratorT && begin, LinesIteratorT && end,
                    float lineWidth, const Color & color, RenderableObject * parent)
                    : GLObject(OpenGLMesh::FromLines(begin, end), 1.0f, lineWidth, 
                    OpenGLShaderSourceDescriptor::DefaultLines, core::Vec3(0, 0, 0), parent),
                    _lines(begin, end), _color(color) {
                    // update color
                    for (auto & v : _mesh.vertices()) {
                        QColor c = MakeQColor(color);
                        v.color4 = core::Vec4f(c.redF(), c.greenF(), c.blueF(), 1.0);
                    }
                }                

            private:
                std::vector<Line3> _lines;
                Color _color;
            };
        }


        // line
        RenderableObject * MakeRenderable(const Line3 & line,
            const DefaultRenderState & state, RenderableObject * parent) {
            return new GLLinesObject(&line, &line + 1, state.lineWidth, state.foregroundColor, parent);
        }

        // lines
        RenderableObject * MakeRenderable(const std::vector<Line3> & lines,
            const DefaultRenderState & state, RenderableObject * parent) {
            return new GLLinesObject(lines.begin(), lines.end(), state.lineWidth, state.foregroundColor, parent);
        }



        namespace {

            OpenGLMesh MeshFromSPPolygons(const std::vector<SpatialProjectedPolygon> & spps, 
                const std::vector<int> & ids){
                OpenGLMesh m;
                for (int id : ids){
                    auto & spp = spps[id];
                    std::vector<Vec3> cs(spp.corners.size());
                    for (int i = 0; i < spp.corners.size(); i++){
                        InfiniteLine3 line(spp.projectionCenter, spp.corners[i] - spp.projectionCenter);
                        cs[i] = IntersectionOfLineAndPlane(line, spp.plane).position;
                    }
                    std::vector<OpenGLMesh::VertHandle> vhandles(cs.size());
                    for (int i = 0; i < cs.size(); i++){
                        vhandles[i] = m.addVertex(Concat(cs[i], 1.0), spp.plane.normal);
                    }
                    m.addPolygon(vhandles);
                }
                return m;
            }

            class GLSpatialProjectedPolygonsObject : public GLObject {
            public:
                explicit GLSpatialProjectedPolygonsObject(const std::vector<SpatialProjectedPolygon> & spps,
                    const std::vector<int> & ids, RenderableObject * parent)
                    : GLObject(MeshFromSPPolygons(spps, ids), 1.0f, 1.0f,
                    OpenGLShaderSourceDescriptor::Panorama, spps[ids.front()].projectionCenter, parent) {
                }
            };

            template <class T, int N>
            struct VecHasher {
                inline size_t operator()(const Vec<T, N> & v) const {
                    size_t h = 0;
                    for (int i = 0; i < N; i++)
                        h = h ^ hasher(v[i]);
                    return h;
                }
                std::hash<T> hasher;
            };

        }


        RenderableObject * MakeRenderable(const std::vector<SpatialProjectedPolygon> & sps,
            const DefaultRenderState & state, RenderableObject * parent){
            std::unordered_map<Vec3, std::vector<int>, VecHasher<double, 3>> idWithSameCenters;
            for (int i = 0; i < sps.size(); i++){
                auto & c = sps[i].projectionCenter;
                idWithSameCenters[c].push_back(i);
            }
            if (idWithSameCenters.size() == 1){
                return new GLSpatialProjectedPolygonsObject(sps, idWithSameCenters.begin()->second, parent);
            }
            else{
                RenderableObject * o = new RenderableObject(parent);
                for (auto & ids : idWithSameCenters){
                    new GLSpatialProjectedPolygonsObject(sps, ids.second, o);
                }
                return o;
            }
        }

        RenderableObject * MakeRenderable(const SpatialProjectedPolygon & sp,
            const DefaultRenderState & state, RenderableObject * parent){
            return MakeRenderable(std::vector<SpatialProjectedPolygon>(1, sp), state, parent);
        }

    }
}