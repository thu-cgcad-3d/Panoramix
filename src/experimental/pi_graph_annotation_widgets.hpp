#pragma once

#include <functional>
#include <QtWidgets>
#include <QtOpenGL>

#include "../gui/scene.hpp"
#include "pi_graph_annotation.hpp"

namespace pano {
    namespace experimental {



        class PIAnnotationWidget : public QGLWidget {
            typedef QGLWidget BaseClass;
        public:
            PIAnnotationWidget(QWidget * parent = nullptr);
            ~PIAnnotationWidget();

            void setCurAnnotation(PIAnnotation * anno);

        protected:
            virtual void paintEvent(QPaintEvent * e) override;
            virtual void mousePressEvent(QMouseEvent * e) override;
            virtual void mouseMoveEvent(QMouseEvent * e) override;
            virtual void mouseReleaseEvent(QMouseEvent * e) override;
            virtual void wheelEvent(QWheelEvent * e) override;
            virtual void keyPressEvent(QKeyEvent * e) override;

        private:
            void clearStroke();
            void acceptAsPolygon(int towardVPId, int alongVPId, bool used);
            void acceptAsOcclusion();

            void rebuildLinesScene();
            void rebuildPolygonScenes();
            void rebuildOcclusionScenes();
            void rebuildStrokeScene();

            void syncToAnnotation();

        private:
            QPoint _lastPos;

            gui::Scene _imageScene;
            gui::Scene _linesScene;
            std::vector<gui::Scene> _polygonScenes;
            std::vector<gui::Scene> _occlusionScenes;
            std::vector<bool> _polygonsDeleted;
            std::vector<bool> _occlusionsDeleted;
            gui::Scene _strokeScene;

            gui::RenderOptions _options;
            PIAnnotation * _anno;

            enum State {
                Idle, CreatingPolygon, CreatingOcclusion
            };
            State _state;
            Chain3 _chain;

            // cur brush
            SegControl _segControl;
        };


      /*  class MainWin : public QMainWindow {
        public:
            MainWin();
            ~MainWin();

            void selectFile(const QString & fname);
            void clear();

            QString imageFileName() const { return _fname; }
            QString annoFileName() const;

        private:
            PIAnnotationWidget * _w;
            PIAnnotation _anno;
            QString _fname;
        };*/


    }
}
