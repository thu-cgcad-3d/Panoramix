#ifndef PANORAMIX_VIS_RENDERABLE_OBJECT_HPP
#define PANORAMIX_VIS_RENDERABLE_OBJECT_HPP

#include <QtOpenGL>
#include "basic_types.hpp"
#include "../core/feature.hpp"

#include "qt_glue.hpp"

#define GL_ALPHA_TEST 0x0BC0

namespace panoramix {
    namespace vis {

        // renderable object
        class RenderableObject : public QObject {
            Q_OBJECT

        public:
            explicit RenderableObject(QObject * parent = nullptr);

        public:
            // render with a given model+view+projection matrix
            virtual void render(RenderModeFlags mode, const QMatrix4x4 & mat) const {}

            // render with given camera and the stored model matrix
            void renderWithCamera(RenderModeFlags mode, const core::PerspectiveCamera & cam) const;

            // bounding box ignoring model matrix
            virtual core::Box3 primaryBoundingBox() const { return core::Box3(); }

            // intersection test ignoring model matrix
            virtual bool intersectsWith(const core::InfiniteLine3 & ray) const { return false; }

            QMatrix4x4 & modelMatrix() { return _modelMat; }
            const QMatrix4x4 & modelMatrix() const { return _modelMat; }

        signals:
            void errorOccored(QString message);

        protected:
            void error(const QString & message);

        protected:
            QMatrix4x4 _modelMat;
        };


        // renderable object list
        using RenderableObjectList = QList<RenderableObject * >;


        // opengl object
        class OpenGLObject : public RenderableObject {
        public:
            explicit OpenGLObject(QObject *parent = 0);
            ~ OpenGLObject();

            void setUpShaders(const OpenGLShaderSource & ss);
            void setUpMesh(const OpenGLMesh & mesh);
            void setUpTexture(const QImage & tex);

            virtual void render(RenderModeFlags mode, const QMatrix4x4 & mat) const override;

        private:
            OpenGLMesh _mesh;
            QOpenGLShaderProgram * _program;
            QOpenGLTexture * _texture;
        };



        // simple lines
        //class OpenGLLines : public QObject, public Renderable {
        //    Q_OBJECT

        //public:
        //    struct Params {
        //        float lineWidth;
        //        Color color;
        //    };

        //    virtual void render(RenderModeFlags mode, const QMatrix4x4 & mat) const;

        //private:
        //    
        //};

    }
}
 
#endif