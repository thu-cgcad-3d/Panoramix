#include "qt_opengl_object.hpp"

namespace panoramix {
    namespace vis {

        void Renderable::render(RenderModeFlags mode, const core::PerspectiveCamera & cam, const QMatrix4x4 & modelMat) const{
            render(mode, MakeQMatrix(cam.viewProjectionMatrix()) * modelMat);
        }

        // opengl object implementation
        OpenGLObject::OpenGLObject(QObject *parent) : QObject(parent), _texture(NULL) {
            _program = new QOpenGLShaderProgram(this);
        }

        OpenGLObject::~OpenGLObject() {
            delete _texture;
        }

        void OpenGLObject::setUpShaders(const OpenGLShaderSource & ss){
            if (!_program->addShaderFromSourceCode(QOpenGLShader::Vertex, 
                QByteArray::fromRawData(ss.vertexShaderSource().data(), ss.vertexShaderSource().size())))
                error(_program->log());
            if (!_program->addShaderFromSourceCode(QOpenGLShader::Fragment, 
                QByteArray::fromRawData(ss.fragmentShaderSource().data(), ss.fragmentShaderSource().size())))
                error(_program->log());
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
            //_texture = new QOpenGLTexture(QImage("F:\\Project.GitHub\\Panoramix\\tests\\data\\panorama\\outdoor\\panohk2.png").mirrored());
            _texture->bind();
            _texture->setMinificationFilter(QOpenGLTexture::Linear);
            _texture->setMagnificationFilter(QOpenGLTexture::Linear);
            _texture->release();
            _program->release();
        }

        void OpenGLObject::setUpMesh(const OpenGLMesh & m){
            _mesh = m;
        }

        void OpenGLObject::render(RenderModeFlags mode, const QMatrix4x4 & mat) const {            

            if (_mesh.vertices().empty())
                return;

            Q_ASSERT(_program->isLinked());
            _program->bind();

            if (_texture && _texture->isCreated())
                _texture->bind(0);

            _program->setUniformValue("matrix", mat);
            _program->setUniformValue("tex", 0);
            _program->setUniformValue("panoramaCenter", QVector3D(0, 0, 0));

            _program->setAttributeArray("position", GL_FLOAT, &(_mesh.vertices().front().position4), 3, sizeof(OpenGLMesh::Vertex));
            _program->setAttributeArray("normal", GL_FLOAT, &(_mesh.vertices().front().normal3), 3, sizeof(OpenGLMesh::Vertex));
            _program->setAttributeArray("color", GL_FLOAT, &(_mesh.vertices().front().color4), 4, sizeof(OpenGLMesh::Vertex));
            _program->setAttributeArray("texCoord", GL_FLOAT, &(_mesh.vertices().front().texCoord2), 2, sizeof(OpenGLMesh::Vertex));

            _program->enableAttributeArray("position");
            _program->enableAttributeArray("normal");
            _program->enableAttributeArray("color");
            _program->enableAttributeArray("texCoord");

            if (mode & RenderModeFlag::Triangles) {
                glDrawElements(GL_TRIANGLES, _mesh.iTriangles().size(), GL_UNSIGNED_INT, _mesh.iTriangles().data());
            } 
            if (mode & RenderModeFlag::Lines) {
                glDrawElements(GL_LINES, _mesh.iLines().size(), GL_UNSIGNED_INT, _mesh.iLines().data());
            }
            if (mode & RenderModeFlag::Points) {
                glDrawElements(GL_POINTS, _mesh.iPoints().size(), GL_UNSIGNED_INT, _mesh.iPoints().data());
            }

            _program->disableAttributeArray("position");
            _program->disableAttributeArray("normal");
            _program->disableAttributeArray("color");
            _program->disableAttributeArray("texCoord");

            _program->release();

        }

        void OpenGLObject::error(const QString & message) {
            qWarning() << message;
            emit errorOccored(message);
        }


    }
}