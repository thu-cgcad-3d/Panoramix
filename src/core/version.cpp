#include "version.hpp"

namespace panoramix {
    namespace core {
 
        std::string Version::toString() const {
            return std::to_string(major) + '.' + std::to_string(minor);
        }

        Version GetVersion() {
            Version v = { PANORAMIX_VERSION_MAJOR, PANORAMIX_VERSION_MINOR };
            return v;
        }
 
    }
}

