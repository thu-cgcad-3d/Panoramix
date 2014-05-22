#pragma once
#ifndef PANORAMIX_DERIV_OPTRAITS_GEOMETRY_HPP
#define PANORAMIX_DERIV_OPTRAITS_GEOMETRY_HPP
 
#include "optraits_base.hpp"

namespace panoramix {
    namespace deriv {

        // map to cross product matrix
        template <class Mat3T, class Vec3T>
        inline Expression<DataStorageType<Mat3T>> makeCross3ProductMatrix(const Expression<Vec3T> & v) {
            using SMat3T = DataStorageType<Mat3T>;
            using SVec3T = DataStorageType<Vec3T>;
            auto mapping = [](ResultType<Vec3T> v) {
                SMat3T mat;
                // a, b: 3x1 vectors
                // a cross b = [a]_{x} * b
                mat << 0, -v(2), v(1),
                    v(2), 0, -v(0),
                    -v(1), v(0), 0;
                return mat;
            };
            auto reverseMapping = [](ResultType<SMat3T> m) {
                return SVec3T(m(2, 1) - m(1, 2), m(0, 2) - m(2, 0), m(1, 0) - m(0, 1));
            };
            return composeMappingFunction<SMat3T, SVec3T>(v.eval(), mapping, reverseMapping, 
                "cross3VecToMat", "cross3MatToVec");
        } 


        // dot product
        namespace {

            template <class T1, class T2>
            struct DotProductTraits 
                : public OpTraitsBase<decltype(specific::DotProduct(std::declval<ResultType<T1>>(), std::declval<ResultType<T2>>())), T1, T2 > {
                inline OutputType value(ResultType<T1> a, ResultType<T2> b) const {
                    return specific::DotProduct(a, b);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T1> a,
                    OriginalAndDerivativeExpression<T2> b) const {
                    a.second = (b.first * sumOfDOutputs).eval();
                    b.second = (a.first * sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "dotProd"; return os; }
            };

        }

        template <class T1, class T2>
        inline auto dotProd(const Expression<T1> & a, const Expression<T2> & b) 
            -> decltype(ComposeExpression(DotProductTraits<T1, T2>(), a, b)) {
            return ComposeExpression(DotProductTraits<T1, T2>(), a, b);
        }

        // norm
        namespace {
            template <class T>
            struct NormTraits : public OpTraitsBase<decltype(std::declval<ResultType<T>>().norm()), T> {
                inline OutputType value(ResultType<T> a) const {
                    return a.norm();
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> a) const {
                    a.second = (a.first * sumOfDOutputs / output).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "norm"; return os; }
            };
        }

        template <class T>
        inline auto norm(const Expression<T> & e)
            ->Expression<decltype(std::declval<ResultType<T>>().norm())> {
            return ComposeExpression(NormTraits<T>(), e);
        }

    }
}

 
#endif