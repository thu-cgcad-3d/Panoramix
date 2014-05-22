#pragma once
#ifndef PANORAMIX_DERIV_OPTRAITS_TRANSPOSE_HPP
#define PANORAMIX_DERIV_OPTRAITS_TRANSPOSE_HPP

#include "optraits_base.hpp"
 
namespace panoramix {
    namespace deriv {

                // transpose
        namespace  {

            template <class T>
            struct TransposeTraits 
                : public OpTraitsBase<decltype(common::Eval(common::GeneralTranspose(std::declval<ResultType<T>>()))), T> {
                inline OutputType value(ResultType<T> from) const {
                    return common::Eval(common::GeneralTranspose(from));
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> from) const {
                    from.second = transpose(sumOfDOutputs).eval().assign<DerivativeType<T>>();
                }
                virtual ostream & toString(ostream & os) const { os << "transpose"; return os; }
            };
        }

        template <class T>
        inline auto transpose(const Expression<T> & e) -> decltype(ComposeExpression(TransposeTraits<T>(), e)) {
            return ComposeExpression(TransposeTraits<T>(), e);
        }


    }
}
 
#endif