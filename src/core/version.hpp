#ifndef PANORAMIX_CORE_VERSION_HPP
#define PANORAMIX_CORE_VERSION_HPP
 
namespace panoramix {
	namespace core {
 
		struct Version {
			const int major, minor;
		};

		// get lib version
		Version GetVersion();
 
	}
}
 
#endif