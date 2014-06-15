#ifndef PANORAMIX_DERIV_TRAITS_HPP
#define PANORAMIX_DERIV_TRAITS_HPP

#include <type_traits>

#include "access.hpp"

namespace panoramix {
    namespace deriv {
        namespace traits {

            typedef std::true_type yes;
            typedef std::false_type no;

            namespace detail {
                //! Used to delay a static_assert until template instantiation
                template <class T>
                struct delay_static_assert : std::false_type {};
            } // namespace detail

            //! Creates a test for whether a non const member function exists
            /*! This creates a class derived from std::integral_constant that will be true if
                the type has the proper member function for the given archive. */
#define DERIV_MAKE_HAS_MEMBER_TEST(name)                                                                                                            \
            namespace detail {                                                                                                                      \
                template <class T, class A>                                                                                                         \
                struct has_member_##name##_impl {                                                                                                   \
                    template <class TT, class AA>                                                                                                   \
                    static auto test(int) -> decltype(panoramix::deriv::access::member_##name(std::declval<AA&>(), std::declval<TT&>()), yes());    \
                    template <class, class>                                                                                                         \
                    static no test(...);                                                                                                            \
                    static const bool value = std::is_same<decltype(test<T, A>(0)), yes>::value;                                                    \
                };                                                                                                                                  \
            } /* end namespace detail */                                                                                                            \
            template <class T, class A>                                                                                                             \
            struct has_member_##name : std::integral_constant<bool, detail::has_member_##name##_impl<T, A>::value>{}

            //! Creates a test for whether a non const non-member function exists
            /*! This creates a class derived from std::integral_constant that will be true if
                the type has the proper non-member function for the given archive. */
#define DERIV_MAKE_HAS_NON_MEMBER_TEST(name)                                                                                  \
            namespace detail {                                                                                                \
                template <class T, class A>                                                                                   \
                struct has_non_member_##name##_impl {                                                                         \
                    template <class TT, class AA>                                                                             \
                    static auto test(int) -> decltype(name(std::declval<AA&>(), std::declval<TT&>()), yes());                 \
                    template <class, class>                                                                                   \
                    static no test(...);                                                                                      \
                    static const bool value = std::is_same<decltype(test<T, A>(0)), yes>::value;                              \
                };                                                                                                            \
            } /* end namespace detail */                                                                                      \
            template <class T, class A>                                                                                       \
            struct has_non_member_##name : std::integral_constant<bool, detail::has_non_member_##name##_impl<T, A>::value>{}



            DERIV_MAKE_HAS_MEMBER_TEST(derEval);
            DERIV_MAKE_HAS_NON_MEMBER_TEST(derEval);
            
            
        }
    }
}
 
#endif