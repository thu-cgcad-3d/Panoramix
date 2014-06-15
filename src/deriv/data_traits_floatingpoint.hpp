#ifndef PANORAMIX_DERIV_DATA_TRAITS_FLOATINGPOINT_HPP
#define PANORAMIX_DERIV_DATA_TRAITS_FLOATINGPOINT_HPP

#include <vector>
#include <numeric>

#include "data_traits.hpp"
 
namespace panoramix {
    namespace deriv {

        DEFINE_DATA_CATEGORY(std::is_floating_point<RemoveAllType<T>>::value, FloatingPointTag);

        template <class T>
        struct DataTraits<T, FloatingPointTag> {
            using StorageType = std::decay_t<T>;
            using ScalarType = std::decay_t<T>;

            static const RoleInProduct roleInProduct = RoleInProduct::Scalar;
            static const bool shouldBeCached = true;

            static StorageType one() { return 1; }
            static StorageType zero() { return 0; }

        };

        namespace common {
            
            // cast
            template <class T1, class T2, 
                SATISFIES(T1, FloatingPointTag), 
                SATISFIES(T2, FloatingPointTag)>
            void Cast(const T1 & from, T2 & to) {
                to = static_cast<T2>(from);
            }

            // eval
            template <class T, SATISFIES(T, FloatingPointTag)>
            std::decay_t<T> Eval(const T & t){
                return t;
            }

            // fill with scalar
            template <class T>
            std::decay_t<T> FillWithScalar(const T & t, const T & s){
                return s;
            }

            // general transpose
            template <class T, SATISFIES(T, FloatingPointTag)>
            std::decay_t<T> GeneralTranspose(const T & t){
                return t;
            }

            // cwise prod
            template <class T1, class T2,
                SATISFIES(T1, FloatingPointTag),
                SATISFIES(T2, FloatingPointTag)>
            auto CWiseProd(const T1 & from, const T2 & to) -> std::decay_t<decltype(from * to)> {
                return from * to;
            }

            // cwise select
            template <class T, class IfT, class ElseT,
                SATISFIES(T, FloatingPointTag),
                SATISFIES(IfT, FloatingPointTag),
                SATISFIES(ElseT, FloatingPointTag)>
            auto CWiseSelect(const T & cond, const IfT & ifval, const ElseT & elseval) -> std::decay_t<decltype(cond > 0 ? ifval : elseval)> {
                return cond > 0 ? ifval : elseval;
            }

            // sum elements
            template <class T, SATISFIES(T, FloatingPointTag)>
            std::decay_t<T> SumElements(const T & d) {
                return d;
            }

            // prod of elements
            template <class T, SATISFIES(T, FloatingPointTag)>
            std::decay_t<T> ProdElements(const T & d) {
                return d;
            }

        }

    }
}
 
#endif