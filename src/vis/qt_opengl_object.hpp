#ifndef PANORAMIX_VIS_QT_OPENGL_OBJECT_HPP
#define PANORAMIX_VIS_QT_OPENGL_OBJECT_HPP

#include <QtOpenGL>
#include "basic_types.hpp"
#include "../core/feature.hpp"

#include "qt_glue.hpp"

#define GL_ALPHA_TEST 0x0BC0

namespace panoramix {
    namespace vis {

        // renderable
        class Renderable {
        public:
            virtual void render(RenderModeFlags mode, const QMatrix4x4 & mat) const = 0;
            void render(RenderModeFlags mode, const core::PerspectiveCamera & cam, const QMatrix4x4 & modelMat) const;
        };

        // opengl object class
        class OpenGLObject : public QObject, public Renderable {
            Q_OBJECT

        public:
            explicit OpenGLObject(QObject *parent = 0);
            ~ OpenGLObject();

            void setUpShaders(const OpenGLShaderSource & ss);
            void setUpMesh(const OpenGLMesh & mesh);
            void setUpTexture(const QImage & tex);

            virtual void render(RenderModeFlags mode, const QMatrix4x4 & mat) const;

        signals:
            void errorOccored(QString message);

        private:
            void error(const QString & message);

        private:
            OpenGLMesh _mesh;
            QOpenGLShaderProgram * _program;
            QOpenGLTexture * _texture;
        };



        // simple lines
        class OpenGLLines : public QObject, public Renderable {
            Q_OBJECT

        public:
            struct Params {
                float lineWidth;
                Color color;
            };

            virtual void render(RenderModeFlags mode, const QMatrix4x4 & mat) const;

        private:
            
        };

    }
}
 
#endif