#ifndef PANORAMIX_DERIV_ACCESS_HPP
#define PANORAMIX_DERIV_ACCESS_HPP
 
namespace panoramix {
    namespace deriv {

        class access {

        	template <class A, class B>
            static auto plus(const A & a, const B & b) -> decltype(a + b) { return a + b; }


        	
        };

    }
}
 
#endif