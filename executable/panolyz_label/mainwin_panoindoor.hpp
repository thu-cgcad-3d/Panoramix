#pragma once

#include <functional>
#include <QtWidgets>
#include <QtOpenGL>

#include "../../src/core/basic_types.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/experimental/rl_graph_annotation.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

namespace panolyz {
    namespace PanoramaIndoor {

        class Widget : public QGLWidget {
            //Q_OBJECT
            typedef QGLWidget BaseClass;
        public:
            Widget(QWidget * parent = nullptr);
            ~Widget();

            void setCurAnnotation(PanoIndoorAnnotation * anno);

        protected:
            virtual void paintEvent(QPaintEvent * e) override;
            virtual void mousePressEvent(QMouseEvent * e) override;
            virtual void mouseMoveEvent(QMouseEvent * e) override;
            virtual void mouseReleaseEvent(QMouseEvent * e) override;
            virtual void wheelEvent(QWheelEvent * e) override;
            virtual void keyPressEvent(QKeyEvent * e) override;

        private:
            void clearStroke();
            void acceptAsPolygon(int towardVPId, int alongVPId, bool clutter);
            void acceptAsOcclusion();

            void rebuildPolygonLineScenes();
            void rebuildStrokeScene();

        private:
            QPoint _lastPos;

            gui::Scene _imageScene;
            std::vector<gui::Scene> _polygonLineScenes;
            gui::Scene _strokeScene;

            gui::RenderOptions _options;
            PanoIndoorAnnotation * _anno;

            enum State {
                Idle, CreatingPolygon, CreatingOcclusion
            };
            State _state;
            Chain3 _chain;

            // cur brush
            int _towardVPId, _alongVPId;
            bool _isClutter;
        };


        class MainWin : public QMainWindow {
            //Q_OBJECT
        public:
            MainWin();
            ~MainWin();

            void selectFile(const QString & fname);
            void clear();

            QString imageFileName() const { return _fname; }
            QString annoFileName() const;

        private:
            Widget * _w;
            PanoIndoorAnnotation _anno;
            QString _fname;
        };


    }
}
