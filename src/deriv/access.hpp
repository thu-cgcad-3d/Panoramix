#ifndef PANORAMIX_DERIV_ACCESS_HPP
#define PANORAMIX_DERIV_ACCESS_HPP
 
namespace panoramix {
    namespace deriv {

        struct access {

        	template <class A, class B>
            static auto plus(const A & a, const B & b) -> decltype(a + b) { return a + b; }

            template <class A, class B>
            static auto minus(const A & a, const B & b) -> decltype(a - b) { return a - b; }

            template <class A, class B>
            static auto cwise_mult(const A & a, const B & b)
                -> std::enable_if_t<std::is_arithmetic<A>::value && 
                std::is_arithmetic<B>::value, decltype(a * b)> { return a * b; }

        	
        };

    }
}
 
#endif