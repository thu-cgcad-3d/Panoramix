#ifndef PANORAMIX_CORE_VERSION_HPP
#define PANORAMIX_CORE_VERSION_HPP

#include <string>
 
namespace panoramix {
    namespace core {
 
        struct Version {
            const int major, minor;
            std::string toString() const;
        };

        // get lib version
        Version GetVersion();
 
    }
}
 
#endif