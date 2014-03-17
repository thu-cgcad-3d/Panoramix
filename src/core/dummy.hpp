#ifndef PANORAMIX_CORE_DUMMY_HPP
#define PANORAMIX_CORE_DUMMY_HPP
 
namespace panoramix {
	namespace core {
 		
 		template <class T>
		struct Dummy {
			T a, b;
		};

		Dummy<int> MakeAnIntDummy();		
 
	}
}
 
#endif