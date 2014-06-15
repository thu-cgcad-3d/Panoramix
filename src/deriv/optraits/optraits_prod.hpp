#pragma once
#ifndef PANORAMIX_DERIV_OPTRAITS_PROD_HPP
#define PANORAMIX_DERIV_OPTRAITS_PROD_HPP
 
#include "optraits_base.hpp"

namespace panoramix {
    namespace deriv {

        // general product
        namespace  {

            // mat * mat
            template <class T1, class T2>
            struct MatrixMatrixProductTraits 
                : public OpTraitsBase<decltype(common::Eval(std::declval<ResultType<T1>>() * std::declval<ResultType<T2>>())), T1, T2> {
                static_assert(DataTraits<T1>::roleInProduct == RoleInProduct::Matrix, "T1 must be a matrix");
                static_assert(DataTraits<T2>::roleInProduct == RoleInProduct::Matrix, "T2 must be a matrix");
                inline OutputType value(ResultType<T1> a, ResultType<T2> b) const {
                    return common::Eval(a * b);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T1> a,
                    OriginalAndDerivativeExpression<T2> b) const {
                    a.second = generalProd(sumOfDOutputs, transpose(b.first)).eval().assign<DerivativeType<T1>>();
                    b.second = generalProd(transpose(a.first), sumOfDOutputs).eval().assign<DerivativeType<T2>>();
                }
                virtual ostream & toString(ostream & os) const { os << "mmProd"; return os; }
            };

            // matrix/array * scalar
            template <class T1, class T2>
            struct MAScalarProductTraits 
                : public OpTraitsBase<decltype(std::declval<ResultType<T1>>() * std::declval<ResultType<T2>>()), T1, T2> {
                static_assert(DataTraits<T1>::roleInProduct == RoleInProduct::Matrix || 
                    DataTraits<T1>::roleInProduct == RoleInProduct::Array, "T1 must be a matrix/array");
                static_assert(DataTraits<T2>::roleInProduct == RoleInProduct::Scalar, "T2 must be a scalar");
                inline OutputType value(ResultType<T1> a, ResultType<T2> b) const {
                    return a * b;
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T1> a,
                    OriginalAndDerivativeExpression<T2> b) const {
                    a.second = generalProd(sumOfDOutputs, b.first).eval();
                    b.second = sumElements(cwiseProd(sumOfDOutputs, a.first)).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "msProd"; return os; }
            };

            // scalar * matrix/array
            template <class T1, class T2>
            struct ScalarMAProductTraits 
                : public OpTraitsBase<decltype(std::declval<ResultType<T1>>() * std::declval<ResultType<T2>>()), T1, T2> {
                static_assert(DataTraits<T1>::roleInProduct == RoleInProduct::Scalar, "T1 must be a scalar");
                static_assert(DataTraits<T1>::roleInProduct == RoleInProduct::Matrix || 
                    DataTraits<T1>::roleInProduct == RoleInProduct::Array, "T2 must be a matrix/array");
                inline OutputType value(ResultType<T1> a, ResultType<T2> b) const {
                    return a * b;
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T1> a,
                    OriginalAndDerivativeExpression<T2> b) const {
                    a.second = sumElements(cwiseProd(sumOfDOutputs, b.first)).eval();
                    b.second = generalProd(sumOfDOutputs, a.first).eval();                    
                }
                virtual ostream & toString(ostream & os) const { os << "smProd"; return os; }
            };

            // array * array
            template <class T1, class T2>
            struct ArrayArrayProductTraits 
                : public OpTraitsBase<decltype(std::declval<ResultType<T1>>() * std::declval<ResultType<T2>>()), T1, T2> {
                static_assert(DataTraits<T1>::roleInProduct == RoleInProduct::Array, "T1 must be an array");
                static_assert(DataTraits<T2>::roleInProduct == RoleInProduct::Array, "T2 must be an array");
                inline OutputType value(ResultType<T1> a, ResultType<T2> b) const {
                    return a * b;
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T1> a,
                    OriginalAndDerivativeExpression<T2> b) const {
                    a.second = cwiseProd(b.first, sumOfDOutputs).eval();
                    b.second = cwiseProd(a.first, sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "aaProd"; return os; }
            };

            // scalar * scalar
            template <class T1, class T2>
            struct ScalarScalarProductTraits
                : public OpTraitsBase<decltype(std::declval<ResultType<T1>>() * std::declval<ResultType<T2>>()), T1, T2> {
                static_assert(DataTraits<T1>::roleInProduct == RoleInProduct::Scalar, "T1 must be a scalar");
                static_assert(DataTraits<T2>::roleInProduct == RoleInProduct::Scalar, "T2 must be a scalar");
                inline OutputType value(ResultType<T1> a, ResultType<T2> b) const {
                    return a * b;
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T1> a,
                    OriginalAndDerivativeExpression<T2> b) const {
                    a.second = cwiseProd(b.first, sumOfDOutputs).eval();
                    b.second = cwiseProd(a.first, sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "ssProd"; return os; }
            };

            template <class T1, class T2, RoleInProduct R1 = DataTraits<T1>::roleInProduct, RoleInProduct R2 = DataTraits<T2>::roleInProduct>
            struct ProductTraits {};

            template <class T1, class T2>
            struct ProductTraits<T1, T2, RoleInProduct::Matrix, RoleInProduct::Matrix> {
                using type = MatrixMatrixProductTraits<T1, T2>;
            };

            template <class T1, class T2>
            struct ProductTraits<T1, T2, RoleInProduct::Matrix, RoleInProduct::Scalar> {
                using type = MAScalarProductTraits<T1, T2>;
            };

            template <class T1, class T2>
            struct ProductTraits<T1, T2, RoleInProduct::Array, RoleInProduct::Scalar> {
                using type = MAScalarProductTraits<T1, T2>;
            };

            template <class T1, class T2>
            struct ProductTraits<T1, T2, RoleInProduct::Scalar, RoleInProduct::Matrix> {
                using type = ScalarMAProductTraits<T1, T2>;
            };

            template <class T1, class T2>
            struct ProductTraits<T1, T2, RoleInProduct::Scalar, RoleInProduct::Array> {
                using type = ScalarMAProductTraits<T1, T2>;
            };

            template <class T1, class T2>
            struct ProductTraits<T1, T2, RoleInProduct::Array, RoleInProduct::Array> {
                using type = ArrayArrayProductTraits<T1, T2>;
            };

            template <class T1, class T2>
            struct ProductTraits<T1, T2, RoleInProduct::Scalar, RoleInProduct::Scalar> {
                using type = ScalarScalarProductTraits<T1, T2>;
            };
        }

        template <class T1, class T2>
        inline auto generalProd(const Expression<T1> & a, const Expression<T2> & b) 
            -> decltype(ComposeExpression(typename ProductTraits<T1, T2>::type(), a, b)) {
            return ComposeExpression(typename ProductTraits<T1, T2>::type(), a, b);
        }

        template <class T1, class T2>
        inline auto operator * (const Expression<T1> & a, const Expression<T2> & b) -> decltype(generalProd(a, b)) {
            return generalProd(a, b);
        }



         // cwise product
        namespace {

            // cwise product traits
            template <class T1, class T2>
            struct CWiseProductTraits 
                : public OpTraitsBase<decltype(common::CWiseProd(std::declval<ResultType<T1>>(), std::declval<ResultType<T2>>())), T1, T2> {
                inline OutputType value(ResultType<T1> a, ResultType<T2> b) const {
                    return common::CWiseProd(a, b);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T1> a,
                    OriginalAndDerivativeExpression<T2> b) const {
                    a.second = cwiseProd(b.first, sumOfDOutputs).eval();
                    b.second = cwiseProd(a.first, sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "cwiseProd"; return os; }
            };

        }

        template <class T1, class T2>
        inline auto cwiseProd(const Expression<T1> & a, const Expression<T2> & b) 
            -> decltype(common::CWiseProd(std::declval<ResultType<T1>>(), std::declval<ResultType<T2>>()), 
                ComposeExpression(CWiseProductTraits<T1, T2>(), a, b)) {
            return ComposeExpression(CWiseProductTraits<T1, T2>(), a, b);
        }



    }
}
 
#endif