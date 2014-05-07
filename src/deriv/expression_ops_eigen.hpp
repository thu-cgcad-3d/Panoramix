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
            static StorageType castFromWithScalarConversion(From from) {
                return from.cast<ScalarType>();
            }

            static StorageType eval(T && t) { return t.eval(); }

            static StorageType fillWithScalar(T && t, ScalarType s) {
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
                using type = decltype(std::declval<ResultType<T>>().array());
            };
            template <class T>
            using MatrixToArrayResultType = typename MatrixToArrayResult<T>::type;

            template <class T>
            struct ArrayToMatrixResult {
                using type = decltype(std::declval<ResultType<T>>().matrix());
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
                    input.second = ArrayToMatrix(sumOfDOutputs).eval();
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
                    input.second = MatrixToArray(sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "array2mat"; return os; }
            };
        }

        template <class MatrixT>
        inline Expression<MatrixToArrayResultType<MatrixT>> MatrixToArray(const Expression<MatrixT> & m) {
            return ComposeExpression(MatrixToArrayTraits<MatrixT>(), m);
        }

        template <class ArrayT>
        inline Expression<ArrayToMatrixResultType<ArrayT>> ArrayToMatrix(const Expression<ArrayT> & a) {
            return ComposeExpression(ArrayToMatrixTraits<ArrayT>(), a);
        }



        // operation of all elements and returns a scalar
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
                    from.second = FillWithScalar(from.first, 1);
                }
                virtual ostream & toString(ostream & os) const { os << "sumElements"; return os; }
            };
        }

        template <class T>
        Expression<DataScalarType<ResultType<T>>> SumElements(const Expression<T> & e){
            return ComposeExpression(SumElementsTraits<T>(), e);
        }


        //// pow
        //namespace  {

        //    using Eigen::pow;

        //    template <class T>
        //    struct PowResult {
        //        using type = decltype(pow(std::declval<ResultType<T>>(), std::declval<DataScalarType<ResultType<T>>>()));
        //    };
        //    template <class T>
        //    using PowResultType = typename PowResult<T>::type;

        //    template <class T>
        //    struct PowTraits : public OpTraitsBase<PowResultType<T>, T> {
        //        inline PowTraits(DataScalarType<T> e) : exponent(e) {}
        //        inline OutputType value(ResultType<T> input) const {
        //            return pow(input, exponent);
        //        }
        //        inline void derivatives(
        //            Expression<OutputType> output,
        //            DerivativeExpression<OutputType> sumOfDOutputs,
        //            OriginalAndDerivativeExpression<T> input) const {
        //            input.second = CWiseProd(pow(input.first, exponent - 1.0) * exponent, sumOfDOutputs).eval();
        //        }
        //        virtual ostream & toString(ostream & os) const { os << "pow[" << exponent << "]"; return os; }
        //        DataScalarType<T> exponent;
        //    };
        //}

        //template <class T>
        //inline Expression<PowResultType<T>> pow(const Expression<T> & a, DataScalarType<T> exponent) {
        //    return ComposeExpression(PowTraits<T>(exponent), a);
        //}

        //template <class T>
        //inline auto operator / (DataScalarType<T> a, const Expression<T> & b) -> decltype(a * pow(b, -1)) {
        //    return a * pow(b, -1.0);
        //}


        //// exp
        //namespace {

        //    using Eigen::exp;

        //    template <class T>
        //    struct ExpResult {
        //        using type = decltype(exp(std::declval<ResultType<T>>()));
        //    };
        //    template <class T>
        //    using ExpResultType = typename ExpResult<T>::type;

        //    template <class T>
        //    struct ExpTraits : public OpTraitsBase<ExpResultType<T>, T> {
        //        inline OutputType value(ResultType<T> input) const {
        //            return exp(input);
        //        }
        //        inline void derivatives(
        //            Expression<OutputType> output,
        //            DerivativeExpression<OutputType> sumOfDOutputs,
        //            OriginalAndDerivativeExpression<T> input) const {
        //            input.second = cwise_product(sumOfDOutputs, output).eval();
        //        }
        //        virtual ostream & toString(ostream & os) const { os << "exp"; return os; }
        //    };

        //}

        //template <class T, IF_EIGEN_ARRAY(T)>
        //inline Expression<ExpResultType<T>> exp(const Expression<T> & a) {
        //    return ComposeExpression(ExpTraits<T>(), a);
        //}


        //// log
        //namespace {

        //    using Eigen::log;

        //    template <class T>
        //    struct LogResult {
        //        using type = decltype(log(std::declval<ResultType<T>>()));
        //    };
        //    template <class T>
        //    using LogResultType = typename LogResult<T>::type;

        //    template <class T>
        //    struct LogTraits : public OpTraitsBase<LogResultType<T>, T> {
        //        inline OutputType value(ResultType<T> input) const {
        //            return log(input);
        //        }
        //        inline void derivatives(
        //            Expression<OutputType> output,
        //            DerivativeExpression<OutputType> sumOfDOutputs,
        //            OriginalAndDerivativeExpression<T> input) const {
        //            input.second = cwise_quotient(sumOfDOutputs, input.first).eval();
        //        }
        //        virtual ostream & toString(ostream & os) const { os << "log"; return os; }
        //    };

        //}

        //template <class T, IF_EIGEN_ARRAY(T)>
        //inline Expression<LogResultType<T>> log(const Expression<T> & a) {
        //    return ComposeExpression(LogTraits<T>(), a);
        //}


        //// sigmoid
        //namespace {

        //    using Eigen::exp;

        //    template <class T>
        //    struct SigmoidResult {
        //        using ScalarType = DataScalarType<T>;
        //        using type = decltype(std::declval<ScalarType>() /
        //            (std::declval<ScalarType>() + exp(-std::declval<ResultType<T>>())));
        //    };
        //    template <class T>
        //    using SigmoidResultType = typename SigmoidResult<T>::type;

        //    template <class T>
        //    struct SigmoidTraits : public OpTraitsBase<SigmoidResultType<T>, T> {
        //        inline OutputType value(ResultType<T> input) const {
        //            return static_cast<DataScalarType<T>>(1) / (static_cast<DataScalarType<T>>(1) + exp(-input));
        //        }
        //        inline void derivatives(
        //            Expression<OutputType> output,
        //            DerivativeExpression<OutputType> sumOfDOutputs,
        //            OriginalAndDerivativeExpression<T> input) const {
        //            input.second = cwise_product(output, 1 - output, sumOfDOutputs).eval();
        //        }
        //        virtual ostream & toString(ostream & os) const { os << "sigmoid"; return os; }
        //    };

        //}

        //template <class T, IF_EIGEN_ARRAY(T)>
        //inline typename SigmoidTraits<T>::OutputExpressionType sigmoid(const Expression<T> & e) {
        //    return ComposeExpression(SigmoidTraits<T>(), e);
        //}



    }
}
 
#endif