#ifndef PANORAMIX_VIS_OPENGL_OBJECT_HPP
#define PANORAMIX_VIS_OPENGL_OBJECT_HPP

#include <QtOpenGL>
#include "../core/basic_types.hpp"
#include "misc.hpp"

#define GL_ALPHA_TEST 0x0BC0

namespace panoramix {
    namespace vis {

        // basic mesh structure for opengl rendering
        struct OpenGLMeshData {
            struct Vertex {
                QVector4D position4;
                QVector3D normal3;
                QVector4D color4;
                QVector2D texCoord2;
                float pointSize1;
                float lineWidth1;
                float a, b;
            };

            using VertHandle = uint32_t;
            using LineHandle = uint32_t;
            using TriangleHandle = uint32_t;

            QVector<Vertex> vertices;
            QVector<VertHandle> iPoints, iLines, iTriangles;

            VertHandle addVertex(const Vertex & v);
            VertHandle addVertex(const QVector4D & p,
                const QVector3D & n = QVector3D(),
                const QVector4D & c = QVector4D(),
                const QVector2D & t = QVector2D());
            LineHandle addLine(VertHandle v1, VertHandle v2);
            LineHandle addIsolatedLine(const Vertex & v1, const Vertex & v2);
            TriangleHandle addTriangle(VertHandle v1, VertHandle v2, VertHandle v3);
            TriangleHandle addIsolatedTriangle(const Vertex & v1, const Vertex & v2, const Vertex & v3);
            void addQuad(VertHandle v1, VertHandle v2, VertHandle v3, VertHandle v4);
            void addPolygon(const QList<VertHandle> & vhs);

            void clear();

            // bounding box
            QPair<QVector3D, QVector3D> boundingBox() const;
        };

        // opengl shader source
        struct OpenGLShaderSource {
            QByteArray vertexShaderSource;
            QByteArray fragmentShaderSource;
        };
        enum class OpenGLShaderSourceName {
            NormalPoints,
            NormalLines,
            NormalTriangles
        };
        OpenGLShaderSource PredefinedShaderSource(OpenGLShaderSourceName name);


        // opengl object class
        class OpenGLObject : public QObject {
            Q_OBJECT

        public:
            explicit OpenGLObject(QObject *parent = 0);
            ~ OpenGLObject();

            void setUpShaders(const OpenGLShaderSource & ss);
            inline void setUpShaders(OpenGLShaderSourceName ssName) { setUpShaders(PredefinedShaderSource(ssName)); }
            void setUpMesh(const OpenGLMeshData & mesh);
            void setUpTexture(const QImage & tex);

            void render(RenderModeFlags mode, const QMatrix4x4 & projection, 
                const QMatrix4x4 & view, const QMatrix4x4 & model);

        signals:
            void errorOccored(QString message);

        private:
            void error(const QString & message);

        private:
            OpenGLMeshData _mesh;
            QOpenGLShaderProgram * _program;
            QOpenGLTexture * _texture;
            GLuint _vertexArrayBuffer;
            GLuint _pointsIndicesBuffer, _linesIndicesBuffer, _trianglesIndicesBuffer;
        };


    }
}
 
#endif