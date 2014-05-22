#pragma once
#ifndef PANORAMIX_DERIV_OPTRAITS_REDUCTION_HPP
#define PANORAMIX_DERIV_OPTRAITS_REDUCTION_HPP

#include "optraits_base.hpp"

namespace panoramix {
    namespace deriv {

         // sum all elements
        namespace  {
            template <class T>
            struct SumElementsTraits : public OpTraitsBase<DataScalarType<ResultType<T>>, T> {
                inline OutputType value(ResultType<T> t) const {
                    return common::SumElements(t);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> from) const {
                    from.second = fillWithScalar(from.first, 1);
                }
                virtual ostream & toString(ostream & os) const { os << "sumElements"; return os; }
            };
        }

        template <class T>
        inline Expression<DataScalarType<ResultType<T>>> sumElements(const Expression<T> & e){
            return ComposeExpression(SumElementsTraits<T>(), e);
        }


    }
}
 
#endif