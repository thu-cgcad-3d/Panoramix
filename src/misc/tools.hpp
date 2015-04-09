#ifndef PANORAMIX_MISC_TOOLS_HPP
#define PANORAMIX_MISC_TOOLS_HPP

#include <algorithm>
#include "../core/basic_types.hpp"

namespace panoramix {
    namespace misc {

        inline std::string Tagify(const std::string & path){
            auto tag = path;
            std::replace(tag.begin(), tag.end(), '.', '_');
            std::replace(tag.begin(), tag.end(), '\\', '.');
            std::replace(tag.begin(), tag.end(), '/', '.');
            std::replace(tag.begin(), tag.end(), ':', '.');
            return tag;
        }

    }
}

#endif