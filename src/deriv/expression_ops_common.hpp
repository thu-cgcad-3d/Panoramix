#ifndef PANORAMIX_DERIV_EXPRESSION_OPS_COMMON_HPP
#define PANORAMIX_DERIV_EXPRESSION_OPS_COMMON_HPP

#include "expression.hpp"
 
namespace panoramix {
    namespace deriv {

        using std::ostream;    


        // fill with scalar
        namespace  {
            template <class T>
            struct FillWithScalarTraits : public OpTraitsBase<DataStorageType<ResultType<T>>, T> {
                using Scalar = DataScalarType<ResultType<T>>;
                inline explicit FillWithScalarTraits(const Scalar & ss) : s(ss) {}
                inline OutputType value(ResultType<T> t) const {
                    return DataTraits<ResultType<T>>::fillWithScalar(std::move(t), s);
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
        inline Expression<DataStorageType<ResultType<T>>> FillWithScalar(Expression<T> t, DataScalarType<ResultType<T>> s) {
            return ComposeExpression(FillWithScalarTraits<T>(s), t);
        }



        // sum
        namespace  {

            template <class T>
            inline T SumAll(T && t){
                return t;
            }

            template <class T, class ... Ts>
            inline auto SumAll(T && t, Ts &&... ts)
                -> decltype(t + SumAll(ts...)) {
                return t + SumAll(ts...);
            }

            template <class ...T>
            struct SumResult;

            template <class T>
            struct SumResult<T> {
                using type = ResultType<T>;
            };

            template <class T, class ... Ts>
            struct SumResult<T, Ts...> : public SumResult<Ts...> {
                using type = decltype(std::declval<ResultType<T>>() + std::declval<typename SumResult<Ts...>::type>());
            };

            template <class ...Ts>
            using SumResultType = typename SumResult<Ts...>::type;

            template <class ...Ts>
            struct SumTraits : public OpTraitsBase<SumResultType<Ts...>, Ts...> {
                inline OutputType value(ResultType<Ts> ... inputs) const {
                    return SumAll(inputs...);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<Ts> ... inputs) const {
                    std::tie(inputs.second ...) =
                        std::tie(
                        sumOfDOutputs.cast<DerivativeType<Ts>>().eval()...);
                }
                virtual ostream & toString(ostream & os) const { os << "sum"; return os; }
            };

        }

        template <class ...Ts>
        inline Expression<SumResultType<Ts...>> GeneralSum(Expression<Ts>... inputs) {
            return ComposeExpression(SumTraits<Ts...>(), inputs...);
        }

        template <class T1, class T2>
        inline auto operator + (Expression<T1> a, Expression<T2> b) -> Expression<decltype(a.result() + b.result())> {
            return ComposeExpression(SumTraits<T1, T2>(), a, b);
        }



        namespace  {

            // homomorphic sums
            template <class T, class EHandleIteratorT>
            inline EHandle HSum2(ExpressionGraph * graph, EHandleIteratorT inpuths) {
                return ComposeExpression(SumTraits<T, T>(),
                    graph->as<T>(*(inpuths)),
                    graph->as<T>(*(inpuths + 1))).eval().handle();
            }
            template <class T, class EHandleIteratorT>
            inline EHandle HSum3(ExpressionGraph * graph, EHandleIteratorT inpuths) {
                return ComposeExpression(SumTraits<T, T, T>(),
                    graph->as<T>(*(inpuths)),
                    graph->as<T>(*(inpuths + 1)),
                    graph->as<T>(*(inpuths + 2))).eval().handle();
            }
            template <class T, class EHandleIteratorT>
            inline EHandle HSum4(ExpressionGraph * graph, EHandleIteratorT inpuths) {
                return ComposeExpression(SumTraits<T, T, T, T>(),
                    graph->as<T>(*(inpuths)),
                    graph->as<T>(*(inpuths + 1)),
                    graph->as<T>(*(inpuths + 2)),
                    graph->as<T>(*(inpuths + 3))).eval().handle();
            }
            template <class T, class EHandleIteratorT>
            inline EHandle HSum5(ExpressionGraph * graph, EHandleIteratorT inpuths) {
                return ComposeExpression(SumTraits<T, T, T, T, T>(),
                    graph->as<T>(*(inpuths)),
                    graph->as<T>(*(inpuths + 1)),
                    graph->as<T>(*(inpuths + 2)),
                    graph->as<T>(*(inpuths + 3)),
                    graph->as<T>(*(inpuths + 4))).eval().handle();
            }

            template <class T>
            inline EHandle HSum(ExpressionGraph * graph, const std::vector<EHandle> & inpuths) {
                std::vector<EHandle> Q(inpuths);
                auto head = Q.begin();
                while (std::distance(head, Q.end()) > 0){
                    EHandle s;
                    switch (std::distance(head, Q.end()))
                    {
                    case 1:
                        return *head;
                    case 2:
                        s = HSum2<T>(graph, head);
                        head += 1;
                        *head = s;
                        break;
                    case 3:
                        s = HSum3<T>(graph, head);
                        head += 2;
                        *head = s;
                        break;
                    case 4:
                        s = HSum4<T>(graph, head);
                        head += 3;
                        *head = s;
                        break;
                    case 5:
                    default:
                        s = HSum5<T>(graph, head);
                        head += 4;
                        *head = s;
                        break;
                    }
                }
                return EHandle();
            }
        }



        // transpose
        namespace  {
            template <class T>
            struct TransposeResult {
                using type = decltype(DataTraits<ResultType<T>>::transpose(std::declval<ResultType<T>>()));
            };
            template <class T>
            using TransposeResultType = typename TransposeResult<T>::type;

            template <class T>
            struct TransposeTraits : public OpTraitsBase<TransposeResultType<T>, T> {
                inline OutputType value(ResultType<T> from) const {
                    return DataTraits<ResultType<T>>::transpose(from);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> from) const {
                    from.second = transpose(sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "transpose"; return os; }
            };
        }

        template <class T>
        Expression<TransposeResultType<T>> transpose(const Expression<T> & e) {
            return ComposeExpression(TransposeTraits<T>(), e);
        }




        // general product
        namespace  {

            // general product
            template <class T>
            inline T ProdAll(T && t){
                return t;
            }

            template <class T, class ... Ts>
            inline auto ProdAll(T && t, Ts &&... ts)
                -> decltype(t * ProdAll(ts...)) {
                return t * ProdAll(ts...);
            }

            template <class ...Ts>
            struct ProductResult;

            template <>
            struct ProductResult<> {};

            template <class T>
            struct ProductResult<T> {
                using type = ResultType<T>;
            };

            template <class T, class ... Ts>
            struct ProductResult<T, Ts...> {
                using type = decltype(std::declval<ResultType<T>>() *
                    std::declval<typename ProductResult<Ts...>::type>());
            };

            template <class ...Ts>
            using ProductResultType = typename ProductResult<Ts...>::type;


            // general product DataTraits
            template <class ...Ts>
            struct ProductTraits : public OpTraitsBase<ProductResultType<Ts...>, Ts...> {
                static_assert(sizeof...(Ts) >= 1, "inputs number must be positive!");
                inline OutputType value(ResultType<Ts> ... inputs) const {
                    return ProdAll(inputs...);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<Ts> ... inputs) const {
                    derivativesUsingSequence(typename SequenceGenerator<InputsNumber>::type(),
                        sumOfDOutputs, inputs...);
                }
                virtual ostream & toString(ostream & os) const { os << "prod"; return os; }

            private:
                template <int ...S>
                inline void derivativesUsingSequence(
                    Sequence<S...> seq,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<Ts> ... inputs) const {
                    auto inputsTuple = std::make_tuple(inputs.first ...);
                    std::tie(inputs.second ...) =
                        std::tie(derivativeWithInputIdx<S>(seq, sumOfDOutputs, inputsTuple)...);
                }

                template <int I, int Idx, class InputExprsTupleT>
                struct _ForEachInputToComposeDInputI {
                    using Type = DerivativeType<TransposeResultType<InputType<I>>>;
                    inline static Expression<Type> expression(
                        DerivativeExpression<OutputType> sumOfDOutputs,
                        const InputExprsTupleT & inputs) {
                        return transpose(std::get<I>(inputs)).eval();
                    }
                };

                template <int Idx, class InputExprsTupleT>
                struct _ForEachInputToComposeDInputI<Idx, Idx, InputExprsTupleT> {
                    using Type = DerivativeType<OutputType>;
                    inline static Expression<Type> expression(
                        DerivativeExpression<OutputType> sumOfDOutputs,
                        const InputExprsTupleT & inputs) {
                        return sumOfDOutputs;
                    }
                };

                template <int InputIdx, class InputExprsTupleT, int ...S>
                inline DerivativeExpression<InputType<InputIdx>> derivativeWithInputIdx(
                    Sequence<S...>,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    const InputExprsTupleT & inputs) const {

                    using TraitsType = ProductTraits<typename
                        _ForEachInputToComposeDInputI<S, InputIdx, InputExprsTupleT>::Type ...>;
                    auto dinput = ComposeExpression(TraitsType(),
                        _ForEachInputToComposeDInputI<S, InputIdx, InputExprsTupleT>::expression(sumOfDOutputs, inputs)...).eval();
                    return ExpressionAssign<DerivativeType<InputType<InputIdx>>>(dinput);
                }
            };

        }

        template <class ...Ts>
        inline Expression<ProductResultType<Ts...>> GeneralProd(const Expression<Ts> & ...inputs) {
            return ComposeExpression(ProductTraits<Ts...>(), inputs...);
        }

        template <class T1, class T2>
        inline auto operator * (const Expression<T1> & a, const Expression<T2> & b) -> decltype(GeneralProd(a, b)) {
            return GeneralProd(a, b);
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
            };
        }

        template <class T>
        inline auto operator - (Expression<T> a) -> Expression<decltype(-a.result())> {
            return ComposeExpression(NegateTraits<T>(), a);
        }


        // subtract
        namespace  {
            template <class T1, class T2>
            struct SubtractTraits : public OpTraitsBase<
                decltype(std::declval<ResultType<T1>>() - std::declval<ResultType<T2>>()), T1, T2> {

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
            };
        }

        template <class T1, class T2>
        inline auto operator - (Expression<T1> a, Expression<T2> b) -> Expression<decltype(a.result() - b.result())> {
            return ComposeExpression(SubtractTraits<T1, T2>(), a, b);
        }


        // mult scalar
        namespace  {
            template <class T>
            struct MultScalar : public OpTraitsBase<
                decltype(std::declval<ResultType<T>>() * std::declval<DataScalarType<ResultType<T>>>()), T> {

                inline explicit MultScalar(DataScalarType<ResultType<T>> ss) : s(ss) {}
                inline OutputType value(ResultType<T> a) const {
                    return a * s;
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> da) const {
                    da.second = (sumOfDOutputs * s).cast<DerivativeType<T>>().eval();
                }

                DataScalarType<ResultType<T>> s;
            };
        }

        template <class T>
        inline auto operator * (Expression<T> a, DataScalarType<ResultType<T>> b)
            ->Expression<decltype(std::declval<ResultType<T>>() * std::declval<DataScalarType<ResultType<T>>>())> {
            return ComposeExpression(MultScalar<T>(b), a);
        }

        template <class T>
        inline auto operator * (DataScalarType<ResultType<T>> b, Expression<T> a)
            ->Expression<decltype(std::declval<ResultType<T>>() * std::declval<DataScalarType<ResultType<T>>>())> {
            return ComposeExpression(MultScalar<T>(b), a);
        }



        // cwise product
        namespace {

            // cwise product
            template <class T>
            inline T CWiseProdAll(T && t){
                return t;
            }

            template <class T, class ... Ts>
            inline auto CWiseProdAll(T && t, Ts && ... ts)
                -> decltype(DataTraits<T>::cwiseProduct(t, CWiseProdAll(ts...))) {
                return DataTraits<T>::cwiseProduct(t, CWiseProdAll(ts...));
            }

            template <class ...Ts>
            struct CWiseProductResult;

            template <>
            struct CWiseProductResult<> {};

            template <class T>
            struct CWiseProductResult<T> {
                using type = ResultType<T>;
            };

            template <class T, class ... Ts>
            struct CWiseProductResult<T, Ts...> {
                using type = decltype(DataTraits<T>::cwiseProduct(std::declval<ResultType<T>>(), 
                    std::declval<typename CWiseProductResult<Ts...>::type>()));
            };

            template <class ...Ts>
            using CWiseProductResultType = typename CWiseProductResult<Ts...>::type;


            // cwise product traits
            template <class ...Ts>
            struct CWiseProductTraits : public OpTraitsBase<CWiseProductResultType<Ts...>, Ts...> {
                static_assert(sizeof...(Ts) >= 1, "inputs number must be positive!");
                inline OutputType value(ResultType<Ts> ... inputs) const {
                    return CWiseProdAll(inputs...);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<Ts> ... inputs) const {
                    derivativesUsingSequence(typename SequenceGenerator<InputsNumber>::type(),
                        sumOfDOutputs, inputs...);
                }
                virtual ostream & toString(ostream & os) const { os << "cwiseProd"; return os; }

            private:
                template <int ...S>
                inline void derivativesUsingSequence(
                    Sequence<S...> seq,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<Ts> ... inputs) const {
                    // tuple of original input expressions
                    auto inputsTuple = std::make_tuple(inputs.first ...);
                    std::tie(inputs.second ...) =
                        std::tie(derivativeWithInputIdx<S>(seq,
                        sumOfDOutputs, inputsTuple)...);
                }

                template <int I, int Idx, class InputExprsTupleT>
                struct _ForEachInputToComposeDInputI {
                    using Type = InputType<I>;
                    inline static Expression<Type> expression(
                        DerivativeExpression<OutputType> sumOfDOutputs,
                        const InputExprsTupleT & inputs) {
                        return std::get<I>(inputs);
                    }
                };

                template <int Idx, class InputExprsTupleT>
                struct _ForEachInputToComposeDInputI<Idx, Idx, InputExprsTupleT> {
                    using Type = DerivativeType<OutputType>;
                    inline static Expression<Type> expression(
                        DerivativeExpression<OutputType> sumOfDOutputs,
                        const InputExprsTupleT & inputs) {
                        return sumOfDOutputs;
                    }
                };

                // compute derivative of each input with split
                template <int InputIdx, class InputExprsTupleT, int ...S>
                inline DerivativeExpression<InputType<InputIdx>> derivativeWithInputIdx(
                    Sequence<S...>,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    const InputExprsTupleT & inputs) const {

                    using TraitsType = CWiseProductTraits<typename _ForEachInputToComposeDInputI<S, InputIdx, InputExprsTupleT>::Type...>;
                    return ComposeExpression(TraitsType(),
                        _ForEachInputToComposeDInputI<S, InputIdx, InputExprsTupleT>::expression(sumOfDOutputs, inputs)...).eval();
                }

            };

        }

        template <class ...Ts>
        inline Expression<CWiseProductResultType<Ts...>> CWiseProd(const Expression<Ts> & ...inputs) {
            return ComposeExpression(CWiseProductTraits<Ts...>(), inputs...);
        }



        // pow
        namespace  {

            using std::pow;
            //using Eigen::pow;

            template <class T>
            struct PowResult {
                using type = decltype(pow(std::declval<ResultType<T>>(), 
                    std::declval<DataScalarType<T>>()));
            };
            template <class T>
            using PowResultType = typename PowResult<T>::type;

            template <class T>
            struct PowTraits : public OpTraitsBase<PowResultType<T>, T> {
                inline PowTraits(DataScalarType<T> e) : exponent(e) {}
                inline OutputType value(ResultType<T> input) const {
                    return pow(input, exponent);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = CWiseProd(pow(input.first, exponent - 1.0) * exponent, sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "pow[" << exponent << "]"; return os; }
                DataScalarType<T> exponent;
            };
        }

        template <class T>
        inline Expression<PowResultType<T>> pow(const Expression<T> & a, DataScalarType<T> exponent) {
            return ComposeExpression(PowTraits<T>(exponent), a);
        }

        template <class T>
        inline auto operator / (DataScalarType<T> a, const Expression<T> & b) -> decltype(a * pow(b, -1)) {
            return a * pow(b, -1.0);
        }

        template <class T1, class T2>
        inline auto CWiseQuotient(const Expression<T1> & a, const Expression<T2> & b) -> decltype(CWiseProd(a, pow(b, -1))) {
            return CWiseProd(a, pow(b, -1));
        }


        // exp
        namespace {
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
                    input.second = CWiseProd(sumOfDOutputs, output).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "exp"; return os; }
            };

        }

        template <class T>
        inline Expression<ExpResultType<T>> exp(const Expression<T> & a) {
            return ComposeExpression(ExpTraits<T>(), a);
        }


        // log
        namespace {            
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
                    input.second = CWiseQuotient(sumOfDOutputs, input.first).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "log"; return os; }
            };

        }

        template <class T>
        inline Expression<LogResultType<T>> log(const Expression<T> & a) {
            return ComposeExpression(LogTraits<T>(), a);
        }


        //// sigmoid
        //namespace {
        //    
        //    using std::exp;
        //    
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

        //template <class T, IF_FLOATING_POINT(T)>
        //inline Expression<SigmoidResultType<T>> sigmoid(const Expression<T> & e) {
        //    return ComposeExpression(SigmoidTraits<T>(), e);
        //}



    }
}
 
#endif