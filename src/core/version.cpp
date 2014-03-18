#include "version.hpp"

namespace panoramix {
	namespace core {
 
		Version GetVersion() {
			Version v = { PANORAMIX_VERSION_MAJOR, PANORAMIX_VERSION_MINOR };
			return v;
		}
 
	}
}

