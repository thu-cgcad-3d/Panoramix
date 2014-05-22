#pragma once
#ifndef PANORAMIX_DERIV_OPTRAITS_SIMPLE_MATH_HPP
#define PANORAMIX_DERIV_OPTRAITS_SIMPLE_MATH_HPP
 
#include "optraits_base.hpp"

namespace panoramix {
    namespace deriv {

         // fill with scalar
        namespace  {
            template <class T>
            struct FillWithScalarTraits : public OpTraitsBase<DataStorageType<T>, T, DataScalarType<T>> {
                using Scalar = DataScalarType<T>;
                inline OutputType value(ResultType<T> t, ResultType<Scalar> s) const {
                    return common::FillWithScalar(t, s);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> t,
                    OriginalAndDerivativeExpression<Scalar> s) const {
                    t.second = DerivativeExpression<T>(); // disconnected
                    s.second = sumOfDOutputs.sum().eval();
                }
                virtual ostream & toString(ostream & os) const { os << "fillWitScalar"; return os; }
            };

            template <class T>
            struct FillWithScalarWhenScalarIsConstTraits : public OpTraitsBase<DataStorageType<T>, T> {
                using Scalar = DataScalarType<T>;
                inline explicit FillWithScalarWhenScalarIsConstTraits(const Scalar & ss) : s(ss) {}
                inline OutputType value(ResultType<T> t) const {
                    return common::FillWithScalar(t, s);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> from) const {
                    from.second = DerivativeExpression<T>(); // disconnected
                }
                virtual ostream & toString(ostream & os) const { os << "fillWitScalar[" << s << "]"; return os; }
                Scalar s;
            };
        }

        template <class T>
        inline Expression<DataStorageType<T>> fillWithScalar(Expression<T> t, Expression<DataScalarType<T>> & s) {
            return ComposeExpression(FillWithScalarTraits<T>(), t, s);
        }

        template <class T>
        inline Expression<DataStorageType<T>> fillWithScalar(Expression<T> t, DataScalarType<T> s) {
            return ComposeExpression(FillWithScalarWhenScalarIsConstTraits<T>(s), t);
        }


        
        // negate
        namespace  {
            template <class T>
            struct NegateTraits : public OpTraitsBase<decltype(-std::declval<ResultType<T>>()), T> {
                inline OutputType value(ResultType<T> a) const {
                    return -a;
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> da) const {
                    da.second = (-sumOfDOutputs).cast<DerivativeType<T>>().eval();
                }
                virtual ostream & toString(ostream & os) const { os << "negate"; return os; }
            };
        }

        template <class T>
        inline auto operator - (Expression<T> a) -> decltype(ComposeExpression(NegateTraits<T>(), a)) {
            return ComposeExpression(NegateTraits<T>(), a);
        }


        // subtract
        namespace  {
            template <class T1, class T2>
            struct SubtractTraits 
                : public OpTraitsBase<decltype(std::declval<ResultType<T1>>() - std::declval<ResultType<T2>>()), T1, T2> {
                inline OutputType value(ResultType<T1> a, ResultType<T2> b) const {
                    return a - b;
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T1> da,
                    OriginalAndDerivativeExpression<T2> db) const {
                    da.second = sumOfDOutputs.cast<DerivativeType<T1>>().eval();
                    db.second = (-sumOfDOutputs).cast<DerivativeType<T2>>().eval();
                }
                virtual ostream & toString(ostream & os) const { os << "subtract"; return os; }
            };
        }

        template <class T1, class T2>
        inline auto operator - (Expression<T1> a, Expression<T2> b) -> decltype(ComposeExpression(SubtractTraits<T1, T2>(), a, b)) {
            return ComposeExpression(SubtractTraits<T1, T2>(), a, b);
        }


        // mult scalar
        namespace  {
            template <class T>
            struct MultScalar 
                : public OpTraitsBase<decltype(std::declval<ResultType<T>>() * std::declval<DataScalarType<ResultType<T>>>()), T> {
                inline explicit MultScalar(DataScalarType<ResultType<T>> ss) : s(ss) {}
                inline OutputType value(ResultType<T> a) const {
                    return a * s;
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> da) const {
                    da.second = (sumOfDOutputs * s).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "multScalar[" << s << "]"; return os; }
                DataScalarType<ResultType<T>> s;
            };
        }

        template <class T>
        inline auto operator * (Expression<T> a, DataScalarType<ResultType<T>> b) -> decltype(ComposeExpression(MultScalar<T>(b), a)) {
            return ComposeExpression(MultScalar<T>(b), a);
        }

        template <class T>
        inline auto operator / (Expression<T> a, DataScalarType<ResultType<T>> b) -> decltype(ComposeExpression(MultScalar<T>(1.0 / b), a)) {
            return ComposeExpression(MultScalar<T>(1.0/b), a);
        }

        template <class T>
        inline auto operator * (DataScalarType<ResultType<T>> b, Expression<T> a) -> decltype(ComposeExpression(MultScalar<T>(b), a)) {
            return ComposeExpression(MultScalar<T>(b), a);
        }

        // plus scalar
        namespace  {
            template <class T>
            struct PlusScalar 
                : public OpTraitsBase<decltype(std::declval<ResultType<T>>() + std::declval<DataScalarType<ResultType<T>>>()), T> {
                inline explicit PlusScalar(DataScalarType<ResultType<T>> ss) : s(ss) {}
                inline OutputType value(ResultType<T> a) const {
                    return a + s;
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> da) const {
                    da.second = (sumOfDOutputs).cast<DerivativeType<T>>().eval();
                }
                virtual ostream & toString(ostream & os) const { os << "plusScalar[" << s << "]"; return os; }
                DataScalarType<ResultType<T>> s;
            };
        }

        template <class T>
        inline auto operator + (Expression<T> a, DataScalarType<ResultType<T>> b) -> decltype(ComposeExpression(PlusScalar<T>(b), a)) {
            return ComposeExpression(PlusScalar<T>(b), a);
        }

        template <class T>
        inline auto operator + (DataScalarType<ResultType<T>> b, Expression<T> a) -> decltype(ComposeExpression(PlusScalar<T>(b), a)) {
            return ComposeExpression(PlusScalar<T>(b), a);
        }

        template <class T>
        inline auto operator - (Expression<T> a, DataScalarType<ResultType<T>> b) -> decltype(ComposeExpression(PlusScalar<T>(-b), a)) {
            return ComposeExpression(PlusScalar<T>(-b), a);
        }

        template <class T>
        inline auto operator - (DataScalarType<ResultType<T>> b, Expression<T> a) -> decltype((-a) + b) {
            return (-a) + b;
        }



       

        // pow
        namespace  {

            using std::pow;

            template <class T>
            struct PowTraits
                : public OpTraitsBase<decltype(pow(std::declval<ResultType<T>>(), std::declval<DataScalarType<T>>())), T> {
                inline PowTraits(DataScalarType<T> e) : exponent(e) {}
                inline OutputType value(ResultType<T> input) const {
                    return pow(input, exponent);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = cwiseProd(pow(input.first, exponent - 1.0) * exponent, sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "pow[" << exponent << "]"; return os; }
                DataScalarType<T> exponent;
            };
        }

        template <class T>
        inline auto pow(const Expression<T> & a, DataScalarType<T> exponent) -> decltype(ComposeExpression(PowTraits<T>(exponent), a)) {
            return ComposeExpression(PowTraits<T>(exponent), a);
        }

        template <class T>
        inline auto operator / (DataScalarType<T> a, const Expression<T> & b) -> decltype(a * pow(b, -1)) {
            return a * pow(b, -1);
        }

        template <class T1, class T2>
        inline auto CWiseQuotient(const Expression<T1> & a, const Expression<T2> & b) -> decltype(cwiseProd(a, pow(b, -1))) {
            return cwiseProd(a, pow(b, -1));
        }

        template <class T1, class T2>
        inline auto operator / (const Expression<T1> & a, const Expression<T2> & b) -> decltype(a * pow(b, -1)) {
            return a * pow(b, -1);
        }


        // exp
        namespace {
            using std::exp;

            template <class T>
            struct ExpTraits : public OpTraitsBase<decltype(exp(std::declval<ResultType<T>>())), T> {
                inline OutputType value(ResultType<T> input) const {
                    return exp(input);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = cwiseProd(sumOfDOutputs, output).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "exp"; return os; }
            };

        }

        template <class T>
        inline auto exp(const Expression<T> & a) -> decltype(ComposeExpression(ExpTraits<T>(), a)){
            return ComposeExpression(ExpTraits<T>(), a);
        }


        // log
        namespace {            
            using std::log;

            template <class T>
            struct LogTraits : public OpTraitsBase<decltype(log(std::declval<ResultType<T>>())), T> {
                inline OutputType value(ResultType<T> input) const {
                    return log(input);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = CWiseQuotient(sumOfDOutputs, input.first).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "log"; return os; }
            };

        }

        template <class T>
        inline auto log(const Expression<T> & a) -> decltype(ComposeExpression(LogTraits<T>(), a)) {
            return ComposeExpression(LogTraits<T>(), a);
        }


        // sigmoid
        namespace {
            template <class T>
            struct SigmoidTraits : public OpTraitsBase<decltype(common::Sigmoid(std::declval<ResultType<T>>())), T> {
                inline OutputType value(ResultType<T> input) const {
                    return common::Sigmoid(input);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = cwiseProd(cwiseProd(output, 1.0 - output), sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "sigmoid"; return os; }
            };
        }

        template <class T>
        inline auto sigmoid(const Expression<T> & e) -> decltype(ComposeExpression(SigmoidTraits<T>(), e)) {
            return ComposeExpression(SigmoidTraits<T>(), e);
        }

        // tanh
        namespace {
            using std::tanh;
            template <class T>
            struct TanhTraits : public OpTraitsBase<decltype(tanh(std::declval<ResultType<T>>())), T> {
                inline OutputType value(ResultType<T> input) const {
                    return tanh(input);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const { // dtanh(x) = (1-tanh(x))^2
                    auto oneMinusOutput = 1.0 - output;
                    input.second = cwiseProd(cwiseProd(oneMinusOutput, oneMinusOutput), sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "tanh"; return os; }
            };
        }

        template <class T>
        inline auto tanh(const Expression<T> & e) -> decltype(ComposeExpression(TanhTraits<T>(), e)) {
            return ComposeExpression(TanhTraits<T>(), e);
        }


        
       




        // abs
        namespace {

            using std::abs;

            template <class T>
            struct AbsTraits : public OpTraitsBase<decltype(abs(std::declval<ResultType<T>>())), T> {
                inline OutputType value(ResultType<T> input) const {
                    return abs(input);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = cwiseSelect(input.first, sumOfDOutputs, -sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "abs"; return os; }
            };
        }

        template <class T>
        inline auto abs(const Expression<T> & e) -> decltype(ComposeExpression(AbsTraits<T>(), e)) {
            return ComposeExpression(AbsTraits<T>(), e);
        }

    }
}
 
#endif