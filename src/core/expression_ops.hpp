#ifndef PANORAMIX_CORE_EXPRESSION_OPS_HPP
#define PANORAMIX_CORE_EXPRESSION_OPS_HPP

#include "expression.hpp"
 
namespace panoramix {
    namespace core { 

        ///////////////////////////////////////////////////////////////////////////////////////
        //// operators
        ///////////////////////////////////////////////////////////////////////////////////////

        // pow
        namespace {

            template <class T>
            struct PowResult {
                using type = decltype(traits<ResultType<T>>::pow(std::declval<ResultType<T>>(),
                std::declval<scalar_type<ResultType<T>>>()));
            };
            template <class T>
            using PowResultType = typename PowResult<T>::type;

            template <class T>
            struct PowTraits : public OpTraitsBase<PowResultType<T>, T> {
                inline PowTraits(scalar_type<T> e) : exponent(e) {}
                inline OutputType value(ResultType<T> input) const {
                    return traits<ResultType<T>>::pow(input, exponent);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = cwise_product(pow(input.first, exponent - 1.0) * exponent, sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "pow[" << exponent << "]"; return os; }
                scalar_type<T> exponent;
            };
        }

        template <class T>
        inline typename PowTraits<T>::OutputExpressionType pow(const Expression<T> & a, scalar_type<T> exponent) {
            return ComposeExpression(PowTraits<T>(exponent), a);
        }

        template <class T>
        inline auto operator / (scalar_type<T> a, const Expression<T> & b) -> decltype(a * pow(b, -1)) {
            return a * pow(b, -1.0);
        }

        template <class T1, class T2>
        inline auto cwise_quotient(const Expression<T1> & a, const Expression<T2> & b)-> decltype(cwise_product(a, pow(b, -1.0))) {
            return cwise_product(a, pow(b, -1.0));
        }
        

        // exp
        namespace {
            
            using Eigen::exp;
            using std::exp;

            template <class T>
            struct ExpResult {
                using type = decltype(exp(std::declval<ResultType<T>>()));
            };
            template <class T>
            using ExpResultType = typename ExpResult<T>::type;

            template <class T>
            struct ExpTraits : public OpTraitsBase<ExpResultType<T>, T> {
                inline OutputType value(ResultType<T> input) const {
                    return exp(input);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = cwise_product(sumOfDOutputs, output).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "exp"; return os; }
            };

        }

        template <class T>
        inline typename ExpTraits<T>::OutputExpressionType exp(const Expression<T> & a) {
            return ComposeExpression(ExpTraits<T>(), a);
        }


        // log
        namespace {

            using Eigen::log;
            using std::log;

            template <class T>
            struct LogResult {
                using type = decltype(log(std::declval<ResultType<T>>()));
            };
            template <class T>
            using LogResultType = typename LogResult<T>::type;

            template <class T>
            struct LogTraits : public OpTraitsBase<LogResultType<T>, T> {
                inline OutputType value(ResultType<T> input) const {
                    return log(input);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = cwise_quotient(sumOfDOutputs, input.first).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "log"; return os; }
            };

        }

        template <class T>
        inline typename LogTraits<T>::OutputExpressionType log(const Expression<T> & a) {
            return ComposeExpression(LogTraits<T>(), a);
        }


        // sigmoid
        namespace {

            using Eigen::exp;
            using std::exp;

            template <class T>
            struct SigmoidResult {
                using ScalarType = scalar_type<T>;
                using type = decltype(std::declval<ScalarType>() / 
                    (std::declval<ScalarType>() + exp(-std::declval<ResultType<T>>())));
            };
            template <class T>
            using SigmoidResultType = typename SigmoidResult<T>::type;

            template <class T>
            struct SigmoidTraits : public OpTraitsBase<SigmoidResultType<T>, T> {
                inline OutputType value(ResultType<T> input) const {
                    return static_cast<scalar_type<T>>(1) / (static_cast<scalar_type<T>>(1) + exp(-input));
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = cwise_product(output, static_cast<scalar_type<T>>(1) - output, sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "sigmoid"; return os; }
            };

        }


    }
}
 
#endif