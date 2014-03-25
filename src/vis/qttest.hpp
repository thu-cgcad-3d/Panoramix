#ifndef PANORAMIX_VIS_QTTEST_HPP
#define PANORAMIX_VIS_QTTEST_HPP

#include <QtWidgets>
#include <QtOpenGL>

#include "../core/feature.hpp"

namespace panoramix {
    namespace vis {

        class TestWidget : public QGLWidget {
        	Q_OBJECT
        public:
        	TestWidget(QWidget * parent = nullptr);
        	~ TestWidget();

        	core::SegmentationExtractor segmenter;
        };

    }
}
 
#endif