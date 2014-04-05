#ifndef PANORAMIX_VIS_QT_RESOURCES_HPP
#define PANORAMIX_VIS_QT_RESOURCES_HPP

#include <memory>

#include "qt_glue.hpp"
#include "qt_opengl_object.hpp"

namespace panoramix {
    namespace vis {

        QApplication* InitGui(int argc, char ** argv);
        QApplication* InitGui();

        void ContinueGui();

    }
}
 
#endif