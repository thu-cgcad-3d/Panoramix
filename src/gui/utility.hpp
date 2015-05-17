#ifndef PANORAMIX_GUI_UTILITY_HPP
#define PANORAMIX_GUI_UTILITY_HPP

#include "basic_types.hpp"

class QWidget;

namespace panoramix {
    namespace gui {

        core::Image PickAnImage(const std::string & dir = std::string(), std::string * picked = nullptr);

        std::vector<core::Image> PickImages(const std::string & dir = std::string());

        std::vector<core::Image> PickAllImagesFromAFolder(const std::string & dir = std::string());

    }
}

#endif