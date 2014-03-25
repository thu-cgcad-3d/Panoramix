#include "visualize3D.hpp"

#include <QtOpenGL>
#include <QtWidgets>

namespace panoramix {
    namespace vis {

        Visualizer3D::Params::Params()
            :
            winName("Image Feature Visualizer"),
            color(255, 255, 255),
            backgroundColor(ColorFromTag(ColorTag::Gray)),
            thickness(1),
            lineType(8),
            colorTableDescriptor(ColorTableDescriptor::AllColors)
        {}



        class VisualizerGUIData : QGLWidget {
        public:
            VisualizerGUIData(QWidget * parent = 0) : QGLWidget(parent){}
        
        protected:
            virtual void initializeGL() {

            }

            virtual void resizeGL(int w, int h) {

            }

            virtual void paintGL() {

            }
            

        public:
            Visualizer3D::Params params;
        };



        Visualizer3D::Visualizer3D(const Params & params, QWidget * parent) 
            : _gui(std::make_unique<VisualizerGUIData>(parent)){
            _gui->params = params;
        }

    }
}