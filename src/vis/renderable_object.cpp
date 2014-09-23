#include <QtOpenGL>

#include "../core/utilities.hpp"

#include "renderable_object.hpp"

#define GL_ALPHA_TEST 0x0BC0

namespace panoramix {
    namespace vis {

        using namespace core;

        RenderableObject::RenderableObject(RenderableObject * parent) : _parent(parent) {
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





        //// opengl object implementation
        //OpenGLObject::OpenGLObject(RenderableObject *parent) : RenderableObject(parent), _texture(NULL) {
        //    _program = new QOpenGLShaderProgram;
        //}

        //OpenGLObject::~OpenGLObject() {
        //    delete _texture;
        //}

        //void OpenGLObject::setUpShaders(const OpenGLShaderSource & ss){
        //    if (!_program->addShaderFromSourceCode(QOpenGLShader::Vertex, 
        //        QByteArray::fromRawData(ss.vertexShaderSource().data(), ss.vertexShaderSource().size())))
        //        error(_program->log());
        //    if (!_program->addShaderFromSourceCode(QOpenGLShader::Fragment, 
        //        QByteArray::fromRawData(ss.fragmentShaderSource().data(), ss.fragmentShaderSource().size())))
        //        error(_program->log());
        //    if (!_program->link() || !_program->bind()){
        //        error(_program->log());
        //        return;
        //    }

        //    Q_ASSERT(_program->isLinked());
        //    qDebug() << _program->log();
        //    _program->release();
        //}


        //void OpenGLObject::setUpTexture(const QImage & tex) {
        //    _program->bind();
        //    Q_ASSERT(_program->isLinked());
        //    if (_texture){
        //        delete _texture;
        //    }
        //    _texture = new QOpenGLTexture(tex.mirrored());
        //    _texture->bind();
        //    _texture->setMinificationFilter(QOpenGLTexture::Linear);
        //    _texture->setMagnificationFilter(QOpenGLTexture::Linear);
        //    _texture->release();
        //    _program->release();
        //}

        //void OpenGLObject::setUpMesh(const OpenGLMesh & m){
        //    _mesh = m;
        //}

        //void OpenGLObject::render(RenderModeFlags mode, const QMatrix4x4 & mat) const {            

        //    if (_mesh.vertices().empty())
        //        return;

        //    Q_ASSERT(_program->isLinked());
        //    _program->bind();

        //    if (_texture && _texture->isCreated())
        //        _texture->bind(0);

        //    _program->setUniformValue("matrix", mat);
        //    _program->setUniformValue("tex", 0);
        //    _program->setUniformValue("panoramaCenter", QVector3D(0, 0, 0));

        //    _program->setAttributeArray("position", GL_FLOAT, &(_mesh.vertices().front().position4), 3, sizeof(OpenGLMesh::Vertex));
        //    _program->setAttributeArray("normal", GL_FLOAT, &(_mesh.vertices().front().normal3), 3, sizeof(OpenGLMesh::Vertex));
        //    _program->setAttributeArray("color", GL_FLOAT, &(_mesh.vertices().front().color4), 4, sizeof(OpenGLMesh::Vertex));
        //    _program->setAttributeArray("texCoord", GL_FLOAT, &(_mesh.vertices().front().texCoord2), 2, sizeof(OpenGLMesh::Vertex));

        //    _program->enableAttributeArray("position");
        //    _program->enableAttributeArray("normal");
        //    _program->enableAttributeArray("color");
        //    _program->enableAttributeArray("texCoord");

        //    if (mode & RenderModeFlag::Triangles) {
        //        glDrawElements(GL_TRIANGLES, _mesh.iTriangles().size(), GL_UNSIGNED_INT, _mesh.iTriangles().data());
        //    } 
        //    if (mode & RenderModeFlag::Lines) {
        //        glDrawElements(GL_LINES, _mesh.iLines().size(), GL_UNSIGNED_INT, _mesh.iLines().data());
        //    }
        //    if (mode & RenderModeFlag::Points) {
        //        glDrawElements(GL_POINTS, _mesh.iPoints().size(), GL_UNSIGNED_INT, _mesh.iPoints().data());
        //    }

        //    _program->disableAttributeArray("position");
        //    _program->disableAttributeArray("normal");
        //    _program->disableAttributeArray("color");
        //    _program->disableAttributeArray("texCoord");

        //    _program->release();

        //}


        //void OpenGLObject::error(const QString & message) {
        //    qWarning() << message;
        //}


        DefaultRenderState::DefaultRenderState() : 
            foregroundColor(0, 0, 0),
            pointSize(10.0f),
            lineWidth(2.0f),
            colorTable(ColorTableDescriptor::AllColors) 
        {}


        namespace {
            class PointsObject : public RenderableObject {
            public:
                template <class PointsIteratorT>
                explicit PointsObject(PointsIteratorT && begin, PointsIteratorT && end,
                    float pointSize, const Color & color, RenderableObject * parent)
                    : RenderableObject(parent), _points(begin, end), _pointSize(pointSize), _color(color) {}

                virtual void render(RenderModeFlags mode, const QMatrix4x4 & mat) const {
                    glPointSize(_pointSize);
                    glBegin(GL_POINTS);
                    for (const Point3 & p : _points) {
                        glColor4dv(_color.val);
                        glVertex3dv(p.val);
                    }
                    glEnd();
                }

                virtual core::Box3 primaryBoundingBox() const { return BoundingBoxOfContainer(_points); }
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
                float _pointSize;
                Color _color;
            };
        }

        // point
        RenderableObject * MakeRenderable(const Point3 & p,
            const DefaultRenderState & param, RenderableObject * parent) {
            return new PointsObject(&p, &p + 1, param.pointSize, param.foregroundColor, parent);
        }

        // points
        RenderableObject * MakeRenderable(const std::vector<Point3> & points,
            const DefaultRenderState & param, RenderableObject * parent) {
            return new PointsObject(points.begin(), points.end(), param.pointSize, param.foregroundColor, parent);
        }




        namespace {
            class LinesObject : public RenderableObject {
            public:
                template <class LinesIteratorT>
                explicit LinesObject(LinesIteratorT && begin, LinesIteratorT && end,
                    float lineWidth, const Color & color, RenderableObject * parent)
                    : RenderableObject(parent), _lines(begin, end), _lineWidth(lineWidth), _color(color) {}

                virtual void render(RenderModeFlags mode, const QMatrix4x4 & mat) const {
                    glLineWidth(_lineWidth);
                    glBegin(GL_LINES);
                    for (const Line3 & l : _lines) {
                        glColor4dv(_color.val);
                        glVertex3dv(l.first.val);
                        glColor4dv(_color.val);
                        glVertex3dv(l.second.val);
                    }
                    glEnd();
                }

                virtual core::Box3 primaryBoundingBox() const { return BoundingBoxOfContainer(_lines); }
                virtual float distanceTo(const core::InfiniteLine3 & ray) const {
                    float distance = std::numeric_limits<float>::max();
                    /*for (auto & l : _lines) {
                        float d = core::DistanceBetweenTwoLines()
                        if (d < distance)
                            distance = d;
                    }*/
                    // TODO
                    return distance;
                }

            private:
                std::vector<Line3> _lines;
                float _lineWidth;
                Color _color;
            };
        }


        // line
        RenderableObject * MakeRenderable(const Line3 & line,
            const DefaultRenderState & state, RenderableObject * parent) {
            return new LinesObject(&line, &line + 1, state.lineWidth, state.foregroundColor, parent);
        }

        // lines
        RenderableObject * MakeRenderable(const std::vector<Line3> & lines,
            const DefaultRenderState & state, RenderableObject * parent) {
            return new LinesObject(lines.begin(), lines.end(), state.lineWidth, state.foregroundColor, parent);
        }




    }
}