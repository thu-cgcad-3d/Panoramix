#pragma once

#include <functional>
#include <QtWidgets>
#include <QtOpenGL>

#include "../gui/scene.hpp"
#include "pi_graph_annotation.hpp"

namespace pano {
    namespace experimental {



        //class PIAnnotationWidget : public QGLWidget {
        //    typedef QGLWidget BaseClass;
        //public:
        //    PIAnnotationWidget(QWidget * parent = nullptr);
        //    ~PIAnnotationWidget();

        //    void setCurAnnotation(PIAnnotation * anno);

        //protected:
        //    virtual void paintEvent(QPaintEvent * e) override;
        //    virtual void mousePressEvent(QMouseEvent * e) override;
        //    virtual void mouseMoveEvent(QMouseEvent * e) override;
        //    virtual void mouseReleaseEvent(QMouseEvent * e) override;
        //    virtual void wheelEvent(QWheelEvent * e) override;
        //    virtual void keyPressEvent(QKeyEvent * e) override;

        //private:
        //    void clearStroke();
        //    void acceptAsPolygon(int towardVPId, int alongVPId, bool used);
        //    void acceptAsOcclusion();
        //    void acceptAsLines();

        //    void rebuildLinesScene();
        //    void rebuildPolygonScenes();
        //    void rebuildOcclusionScenes();
        //    void rebuildStrokeScene();

        //private:
        //    QPoint _lastPos;

        //    gui::Scene _imageScene;
        //    gui::Scene _linesScene;
        //    std::vector<gui::Scene> _polygonScenes;
        //    std::vector<gui::Scene> _occlusionScenes;
        //    std::vector<bool> _polygonsDeleted;
        //    std::vector<bool> _occlusionsDeleted;
        //    gui::Scene _strokeScene;

        //    gui::RenderOptions _options;
        //    PIAnnotation * _anno;

        //    enum State {
        //        Idle, CreatingPolygon, CreatingOcclusion, CreatingLine
        //    };
        //    State _state;
        //    Chain3 _chain;

        //    // cur brush
        //    SegControl _segControl;

        //    bool _showPolygons;
        //    bool _showLines;
        //    bool _showOcclusions;
        //    bool _showVPs;
        //};





        class PILayoutAnnotationWidget : public QGLWidget {
        public:
            PILayoutAnnotationWidget(QWidget * parent = nullptr);
            ~PILayoutAnnotationWidget();

            void setCurAnnotation(PILayoutAnnotation * anno);

        protected:
            virtual void paintEvent(QPaintEvent * e) override;
            virtual void mousePressEvent(QMouseEvent * e) override;
            virtual void mouseMoveEvent(QMouseEvent * e) override;
            virtual void mouseReleaseEvent(QMouseEvent * e) override;
            virtual void wheelEvent(QWheelEvent * e) override;
            virtual void keyPressEvent(QKeyEvent * e) override;

        private:
            void clearStroke();

            void rebuildLayoutScene();
            void rebuildCluttersScene();
            void rebuildStrokeScene();

            void acceptClutter();

        private:
            QPoint _lastPos;

            gui::Scene _imageScene;

            gui::Scene _layoutScene;
            std::vector<Decorated<gui::Colored<Point3>, int>> _cornerPoints;
            std::vector<Decorated<gui::Colored<Line3>, int>> _borderLines;
            std::vector<Decorated<gui::Colored<Polygon3>, int>> _facePolygons;
            
            gui::Scene _cluttersScene;
            std::vector<Decorated<Polygon3, int>> _clutterPolygons;

            gui::Scene _vpsScene;
            gui::Scene _strokeScene;

            gui::RenderOptions _options;
            PILayoutAnnotation * _anno;

            enum State {
                Idle, DrawingBorder, DrawingClutter
            };
            State _state;
            int _lastHitCornerId;
            Chain3 _stroke;

            int _cornerClicked, _borderClicked, _faceClicked;

            bool _showLayouts;
            bool _showClutters;
            bool _showVPs;
        };


    }
}
