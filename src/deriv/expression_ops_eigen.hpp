#ifndef PANORAMIX_CORE_EXPRESSION_OPS_EIGEN_HPP
#define PANORAMIX_CORE_EXPRESSION_OPS_EIGEN_HPP

#include <Eigen/Core>
#include <Eigen/Dense>

#include "expression.hpp"
#include "expression_ops_common.hpp"
 
namespace panoramix {
    namespace deriv {

        // eigen classes
        namespace {

            template <class T>
            struct is_eigen_dense
                : public std::is_base_of<Eigen::EigenBase<T>, T>
            {};

            template <class Derived>
            struct is_eigen_dense<Eigen::DenseBase<Derived>> : public std::true_type{};

            template <class T>
            struct is_eigen_matrix
                : public std::is_base_of<Eigen::MatrixBase<T>, T>
            {};

            template <class Derived>
            struct is_eigen_matrix<Eigen::MatrixBase<Derived>> : public std::true_type{};

            template <class T>
            struct is_eigen_array
                : public std::is_base_of<Eigen::ArrayBase<T>, T>
            {};

            template <class Derived>
            struct is_eigen_array<Eigen::ArrayBase<Derived>> : public std::true_type{};

        }

        DEFINE_DATA_CATEGORY(is_eigen_dense<RemoveAllType<T>>::value &&
            !is_eigen_array<RemoveAllType<T>>::value &&
            !is_eigen_matrix<RemoveAllType<T>>::value,
            EigenDenseTag);
        DEFINE_DATA_CATEGORY(is_eigen_array<RemoveAllType<T>>::value, EigenArrayTag);
        DEFINE_DATA_CATEGORY(is_eigen_matrix<RemoveAllType<T>>::value, EigenMatrixTag);

        template <class T>
        struct DataTraits<T, EigenDenseTag> {

            using StorageType = std::decay_t<decltype(std::declval<T>().eval())>;
            using ScalarType = typename StorageType::Scalar;

            static const bool shouldBeCached = std::is_same<T, StorageType>::value;

            template <class From>
            static StorageType castFromWithScalarConversion(const From & from) {
                return from.cast<ScalarType>();
            }

            static StorageType eval(const T & t) { 
                return t.eval(); 
            }

            static StorageType fillWithScalar(const T & t, ScalarType s) {
                StorageType tt = t;
                tt.setConstant(s);
                return tt;
            }

            static StorageType one() { 
                assert(false && "cannot get one from eigen dense (which is not a scalar) type!"); 
                throw ""; 
            }

            static StorageType zero() {
                assert(false && "cannot get zero from eigen dense (which is not a scalar) type!");
                throw ""; 
            }

        };


        template <class T>
        struct DataTraits<T, EigenArrayTag> : public DataTraits<T, EigenDenseTag>{

            template <class K>
            static auto cwiseProduct(const T & a, K && b) -> decltype(a.cwiseProduct(b))  {
                return a.cwiseProduct(b);
            }
            
            static T transpose(const T & t) {
                return t;
            }

        };

        template <class T>
        struct DataTraits<T, EigenMatrixTag> : public DataTraits<T, EigenDenseTag>{

            template <class K>
            static auto cwiseProduct(const T & a, K && b) -> decltype(a.cwiseProduct(b))  {
                return a.cwiseProduct(b);
            }

            static auto transpose(const T & t) -> decltype(t.transpose()) { 
                return t.transpose(); 
            }

        };




        // matrix to array
        // array to matrix
        namespace  {

            template <class T>
            struct MatrixToArrayResult {
                using type = DataStorageType<decltype(std::declval<ResultType<T>>().array())>;
            };
            template <class T>
            using MatrixToArrayResultType = typename MatrixToArrayResult<T>::type;

            template <class T>
            struct ArrayToMatrixResult {
                using type = DataStorageType<decltype(std::declval<ResultType<T>>().matrix())>;
            };
            template <class T>
            using ArrayToMatrixResultType = typename ArrayToMatrixResult<T>::type;

            // matrix to array
            template <class T>
            struct MatrixToArrayTraits : public OpTraitsBase<MatrixToArrayResultType<T>, T> {
                inline OutputType value(ResultType<T> input) const {
                    return input.array();
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = arrayToMatrix(sumOfDOutputs);
                }
                virtual ostream & toString(ostream & os) const { os << "mat2array"; return os; }
            };

            // array to matrix
            template <class T>
            struct ArrayToMatrixTraits : public OpTraitsBase<ArrayToMatrixResultType<T>, T> {
                inline OutputType value(ResultType<T> input) const {
                    return input.matrix();
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = matrixToArray(sumOfDOutputs);
                }
                virtual ostream & toString(ostream & os) const { os << "array2mat"; return os; }
            };
        }

        template <class MatrixT>
        inline Expression<MatrixToArrayResultType<MatrixT>> matrixToArray(const Expression<MatrixT> & m) {
            return ComposeExpression(MatrixToArrayTraits<MatrixT>(), m);
        }

        template <class ArrayT>
        inline Expression<ArrayToMatrixResultType<ArrayT>> arrayToMatrix(const Expression<ArrayT> & a) {
            return ComposeExpression(ArrayToMatrixTraits<ArrayT>(), a);
        }



        // sum all elements
        namespace  {
            template <class T>
            struct SumElementsTraits : public OpTraitsBase<DataScalarType<ResultType<T>>, T> {
                inline OutputType value(ResultType<T> t) const {
                    return t.sum();
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


        // inverse matrix
        namespace {
            template <class MatrixT>
            struct InverseMatrixTraits : public OpTraitsBase<DataStorageType<MatrixT>, MatrixT> {
                static_assert(is_eigen_matrix <RemoveAllType<MatrixT>>::value, "MatrixT MUST be an eigen matrix!");
                inline OutputType value(ResultType<MatrixT> t) const {
                    return t.inverse();
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