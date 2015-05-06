#ifndef PANORAMIX_GUI_UTILITIES_HPP
#define PANORAMIX_GUI_UTILITIES_HPP

#include "basic_types.hpp"

class QWidget;

namespace panoramix {
    namespace gui {

        core::Image PickAnImage(const std::string & dir = std::string());

        std::vector<core::Image> PickImages(const std::string & dir = std::string());



    }
}

#endif