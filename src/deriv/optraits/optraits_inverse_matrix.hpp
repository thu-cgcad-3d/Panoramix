#pragma once
#ifndef PANORAMIX_DERIV_OPTRAITS_INVERSE_MATRIX_HPP
#define PANORAMIX_DERIV_OPTRAITS_INVERSE_MATRIX_HPP

#include "optraits_base.hpp"

namespace panoramix {
    namespace deriv {

        // inverse matrix
        namespace {
            template <class MatrixT>
            struct InverseMatrixTraits : public OpTraitsBase<DataStorageType<MatrixT>, MatrixT> {
                inline OutputType value(ResultType<MatrixT> t) const {
                    return specific::InverseMatrix(t);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<MatrixT> from) const {
                    from.second = transpose(-output * transpose(sumOfDOutputs) * output)
                        .eval().assign<DerivativeType<MatrixT>>();
                }
                virtual ostream & toString(ostream & os) const { os << "inverseMatrix"; return os; }
            };
        }

        template <class T>
        inline Expression<DataStorageType<T>> inverseMatrix(const Expression<T> & m){
            return ComposeExpression(InverseMatrixTraits<T>(), m);
        }

    }
}
 
#endif