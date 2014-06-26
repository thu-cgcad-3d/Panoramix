#pragma once
#ifndef PANORAMIX_DERIV_OPTRAITS_EIGEN_HPP
#define PANORAMIX_DERIV_OPTRAITS_EIGEN_HPP

#include "optraits_base.hpp"

namespace panoramix {
    namespace deriv {

        // matrix to array
        // array to matrix
        namespace  {

            // matrix to array
            template <class T>
            struct MatrixToArrayTraits 
                : public OpTraitsBase<decltype(common::Eval(std::declval<ResultType<T>>().array())), T> {
                inline OutputType value(ResultType<T> input) const {
                    return input.array();
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = arrayToMatrix(sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "matrixToArray"; return os; }
            };

            // array to matrix
            template <class T>
            struct ArrayToMatrixTraits 
                : public OpTraitsBase<decltype(common::Eval(std::declval<ResultType<T>>().matrix())), T> {
                inline OutputType value(ResultType<T> input) const {
                    return input.matrix();
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = matrixToArray(sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "arrayToMatrix"; return os; }
            };
        }

        template <class MatrixT>
        inline auto matrixToArray(const Expression<MatrixT> & m) 
            -> decltype(ComposeExpression(MatrixToArrayTraits<MatrixT>(), m)) {
            return ComposeExpression(MatrixToArrayTraits<MatrixT>(), m);
        }

        template <class ArrayT>
        inline auto arrayToMatrix(const Expression<ArrayT> & a) 
            -> decltype(ComposeExpression(ArrayToMatrixTraits<ArrayT>(), a)) {
            return ComposeExpression(ArrayToMatrixTraits<ArrayT>(), a);
        }



        // replicate
        namespace {

            template <class Derived>
            inline auto Replicate(const Eigen::DenseBase<Derived> & d, 
                typename Eigen::DenseBase<Derived>::Index rowFactor, typename Eigen::DenseBase<Derived>::Index colFactor)
                -> decltype(d.replicate(rowFactor, colFactor)) {
                return d.replicate(rowFactor, colFactor);
            }

            template <class Derived>
            inline auto FoldSum(const Eigen::DenseBase<Derived> & d,
                typename Eigen::DenseBase<Derived>::Index rowFactor, typename Eigen::DenseBase<Derived>::Index colFactor)
                -> decltype(d.block(0, 0, 0, 0).eval()) {
                assert(d.rows() % rowFactor == 0 && "row num incompitable with rowFactor!");
                assert(d.cols() % colFactor == 0 && "col num incompitable with colFactor!");
                typename Eigen::DenseBase<Derived>::Index left = 0, top = 0, rows = d.rows() / rowFactor, cols = d.cols() / colFactor;
                auto result = d.block(0, 0, rows, cols).eval();
                for (typename Eigen::DenseBase<Derived>::Index i = 0; i < d.rows(); i += rows){
                    for (typename Eigen::DenseBase<Derived>::Index j = 0; j < d.cols(); j += cols){
                        if (i == 0 && j == 0) continue;
                        result += d.block(i, j, rows, cols);
                    }
                }
                return result;
            }


            template <class T>
            struct ReplicateTraits 
                : public OpTraitsBase<decltype(Replicate(std::declval<ResultType<T>>(), 0, 0)), T> {
                inline explicit ReplicateTraits(int r, int c) : rowFactor(r), colFactor(c) {}
                inline OutputType value(ResultType<T> input) const {
                    return Replicate(input, rowFactor, colFactor);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = foldSum(sumOfDOutputs, rowFactor, colFactor).eval().assign<DerivativeType<T>>();
                }
                virtual ostream & toString(ostream & os) const { os << "replicate"; return os; }
                int rowFactor, colFactor;
            };

            template <class T, class RowPred, class ColPred>
            struct replicateTraits
                : public OpTraitsBase<decltype(Replicate(std::declval<ResultType<T>>(), std::declval<RowPred>()(), std::declval<ColPred>()())), T> {
                inline explicit replicateTraits(RowPred r, ColPred c) : rowFactor(r), colFactor(c) {}
                inline OutputType value(ResultType<T> input) const {
                    return Replicate(input, rowFactor(), colFactor());
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = foldSumWithPred(sumOfDOutputs, rowFactor, colFactor).eval().assign<DerivativeType<T>>();
                }
                virtual ostream & toString(ostream & os) const { os << "replicate"; return os; }
                RowPred rowFactor;
                ColPred colFactor;
            };
            
            template <class T>
            struct FoldSumTraits
                : public OpTraitsBase<decltype(FoldSum(std::declval<ResultType<T>>(), 0, 0)), T> {
                inline explicit FoldSumTraits(int r, int c) : rowFactor(r), colFactor(c) {}
                inline OutputType value(ResultType<T> input) const {
                    return FoldSum(input, rowFactor, colFactor);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = replicate(sumOfDOutputs, rowFactor, colFactor).eval().assign<DerivativeType<T>>();
                }
                virtual ostream & toString(ostream & os) const { os << "disreplicate"; return os; }
                int rowFactor, colFactor;
            };

            template <class T, class RowPred, class ColPred>
            struct FoldSumWithPredTraits
                : public OpTraitsBase<decltype(FoldSum(std::declval<ResultType<T>>(), std::declval<RowPred>()(), std::declval<ColPred>()())), T> {
                inline explicit FoldSumWithPredTraits(RowPred r, ColPred c) : rowFactor(r), colFactor(c) {}
                inline OutputType value(ResultType<T> input) const {
                    return FoldSum(input, rowFactor(), colFactor());
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = replicate(sumOfDOutputs, rowFactor, colFactor).eval().assign<DerivativeType<T>>();
                }
                virtual ostream & toString(ostream & os) const { os << "disreplicate"; return os; }
                RowPred rowFactor;
                ColPred colFactor;
            };

        }


        template <class T>
        inline auto replicate(const Expression<T> & a, int rowFactor, int colFactor)
            -> decltype(Replicate(std::declval<T>(), rowFactor, colFactor), ComposeExpression(ReplicateTraits<T>(rowFactor, colFactor), a)) {
            return ComposeExpression(ReplicateTraits<T>(rowFactor, colFactor), a);
        }

        template <class T, class RowPred, class ColPred>
        inline auto replicate(const Expression<T> & a, RowPred rowFactor, ColPred colFactor)
            -> decltype(Replicate(std::declval<T>(), rowFactor(), colFactor()), ComposeExpression(replicateTraits<T, RowPred, ColPred>(rowFactor, colFactor), a)) {
            return ComposeExpression(replicateTraits<T, RowPred, ColPred>(rowFactor, colFactor), a);
        }


        template <class T>
        inline auto foldSum(const Expression<T> & a, int rowFactor, int colFactor)
            -> decltype(FoldSum(std::declval<T>(), rowFactor, colFactor), ComposeExpression(FoldSumTraits<T>(rowFactor, colFactor), a)) {
            return ComposeExpression(FoldSumTraits<T>(rowFactor, colFactor), a);
        }

        template <class T, class RowPred, class ColPred>
        inline auto foldSumWithPred(const Expression<T> & a, RowPred rowFactor, ColPred colFactor)
            -> decltype(FoldSum(std::declval<T>(), rowFactor(), colFactor()), ComposeExpression(FoldSumWithPredTraits<T, RowPred, ColPred>(rowFactor, colFactor), a)) {
            return ComposeExpression(FoldSumWithPredTraits<T, RowPred, ColPred>(rowFactor, colFactor), a);
        }

    }
}
 
#endif