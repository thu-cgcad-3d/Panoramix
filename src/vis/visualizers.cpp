#include <QtOpenGL>

#include "qt_glue.hpp"
#include "visualizers.hpp"

#define GL_ALPHA_TEST 0x0BC0

namespace panoramix {
    namespace vis {

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



        ResourcePtr MakeTextureResource(const core::Image & image){
            struct TextureResource : Resource {
                inline TextureResource(const core::Image & im) : 
                    image(MakeQImage(im)),
                    texture(new QOpenGLTexture(QOpenGLTexture::Target2D)) {}
                virtual bool isNull() const override { return image.isNull(); }
                virtual bool initialize() override {
                    if (!texture->isCreated()){
                        texture->create();
                    }
                    if (!texture->isCreated()){
                        return false;
                    }
                    Q_ASSERT(texture->textureId());
                    texture->setData(image.mirrored());
                    texture->setMinificationFilter(QOpenGLTexture::Linear);
                    texture->setMagnificationFilter(QOpenGLTexture::Linear);
                    texture->release();
                    return true;
                }
                virtual bool bind() override { texture->bind(0); return texture->isBound(); }
                virtual ~TextureResource() { delete texture; }
                QImage image;
                QOpenGLTexture * texture;
            };

            return std::make_unique<TextureResource>(image);
        }



        struct VisualObjectInternal {
            QOpenGLShaderProgram * program;            
            inline VisualObjectInternal() :
                program(new QOpenGLShaderProgram) {}
            inline ~VisualObjectInternal(){
                delete program;
            }
        };

        VisualObject::VisualObject(const OpenGLShaderSource & shaderSource, VisualObject * parent)
            : _shaderSource(shaderSource), _parent(parent), _modelMat(core::Mat4::eye()) {
            _internal = new VisualObjectInternal;
            if (parent){
                parent->_children.push_back(this);
            }
        }

        VisualObject::~VisualObject(){
            delete _internal;
        }

        void VisualObject::initialize() const {
            auto vo = static_cast<VisualObjectInternal*>(_internal);
            auto program = vo->program;

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
                if (!res->initialize())
                    qDebug() << (program->log());
            }           

            program->release();
        }

        void VisualObject::render(const Options & options, const core::Mat4 & mat) const {
            if (_mesh.vertices.empty())
                return;

            auto vo = static_cast<VisualObjectInternal*>(_internal);
            auto program = vo->program;

            Q_ASSERT(program->isLinked());
            program->bind();

            // bind resources
            for (auto & res : _resources){
                if (!res->bind())
                    qDebug() << (program->log());
            }            

            glEnable(GL_PROGRAM_POINT_SIZE);
            glEnable(GL_POINT_SPRITE);

            program->setUniformValue("useUniformLineWidth", options.lineWidth.enabled);
            if (options.lineWidth.enabled){
                glLineWidth(options.lineWidth.component);
            }
            program->setUniformValue("useUniformPointSize", options.pointSize.enabled);
            if (options.pointSize.enabled){
                glPointSize(options.pointSize.component);
            }

            if (options.projectionCenter.enabled)
                program->setUniformValue("panoramaCenter", MakeQVec(options.projectionCenter.component));

            program->setUniformValue("matrix", MakeQMatrix(mat));
            program->setUniformValue("tex", 0);            

            SetAttributeArrayWithTriMeshVertices(program, "position", _mesh.vertices.front().position4);
            SetAttributeArrayWithTriMeshVertices(program, "normal", _mesh.vertices.front().normal3);
            SetAttributeArrayWithTriMeshVertices(program, "color", _mesh.vertices.front().color4);
            SetAttributeArrayWithTriMeshVertices(program, "texCoord", _mesh.vertices.front().texCoord2);
            SetAttributeArrayWithTriMeshVertices(program, "pointSize", _mesh.vertices.front().pointSize1);

            program->enableAttributeArray("position");
            program->enableAttributeArray("normal");
            program->enableAttributeArray("color");
            program->enableAttributeArray("texCoord");
            program->enableAttributeArray("pointSize");

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
            program->disableAttributeArray("pointSize");
            
            program->release();

        }

        void VisualObject::renderWithCamera(const Options & options, const core::PerspectiveCamera & cam) const {
            render(options, cam.viewProjectionMatrix() * _modelMat);
        }





    }
}