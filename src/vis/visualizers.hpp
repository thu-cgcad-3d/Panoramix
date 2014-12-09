#ifndef PANORAMIX_VIS_VISUALIZERS_HPP
#define PANORAMIX_VIS_VISUALIZERS_HPP

#include "../core/feature.hpp"
#include "../core/cameras.hpp"
#include "basic_types.hpp"

namespace panoramix {
    namespace vis {

        struct Options {
            std::string winName;
            Color backgroundColor;
            RenderModeFlags renderMode;
        };

        

    }
}
 
#endif