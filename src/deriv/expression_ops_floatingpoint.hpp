#ifndef PANORAMIX_CORE_EXPRESSION_OPS_FLOATINGPOINT_HPP
#define PANORAMIX_CORE_EXPRESSION_OPS_FLOATINGPOINT_HPP

#include <vector>
#include <numeric>

#include "expression.hpp"
#include "expression_ops_common.hpp"
 
namespace panoramix {
    namespace deriv {    

        DEFINE_DATA_CATEGORY(std::is_floating_point<RemoveAllType<T>>::value, FloatingPointTag);

        template <class T>
        struct DataTraits<T, FloatingPointTag> {
            using StorageType = std::decay_t<T>;
            using ScalarType = std::decay_t<T>;

            static const bool shouldBeCached = true;

            template <class From>
            static StorageType castFromWithScalarConversion(From from) {
                return static_cast<StorageType>(from);
            }

            static StorageType eval(T t) { return t; }

            static StorageType fillWithScalar(T t, ScalarType s) { return s; }

            static StorageType one() { return 1; }

            static StorageType zero() { return 0; }

            static StorageType transpose(T t) { return t; }

            template <class K>
            static auto cwiseProduct(T a, K b) -> decltype(a * b) { return a * b; }
        };

    }
}
 
#endif