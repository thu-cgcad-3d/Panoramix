#include "cache.hpp"

namespace pano {
    namespace misc {

        std::string Tagify(const std::string & path) {
            auto tag = path;
            std::replace(tag.begin(), tag.end(), '.', '_');
            std::replace(tag.begin(), tag.end(), '\\', '.');
            std::replace(tag.begin(), tag.end(), '/', '.');
            std::replace(tag.begin(), tag.end(), ':', '.');
            return tag;
        }

    }
}