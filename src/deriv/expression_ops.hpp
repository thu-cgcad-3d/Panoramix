#ifndef PANORAMIX_DERIV_EXPRESSION_OPS_HPP
#define PANORAMIX_DERIV_EXPRESSION_OPS_HPP

#include "expression.hpp"
 
namespace panoramix {
    namespace deriv {

        using std::ostream;



        ////////////////////////////////////////////////////////////////////////////
        // the functional op wrapper
        template <class OutputT, class ...InputTs>
        struct OpTraitsBase {
            using OutputType = OutputT;
            using OutputExpressionType = Expression<OutputT>;

            static const int InputsNumber = sizeof...(InputTs);
            using InputIndices = typename SequenceGenerator<sizeof...(InputTs)>::type;
            using InputTuple = std::tuple<InputTs...>;

            template <int I>
            struct InputTypeStruct {
                using type = typename std::tuple_element<I, InputTuple>::type;
            };
            template <int I>
            using InputType = typename InputTypeStruct<I>::type;

            // to string
            virtual ostream & toString(ostream & os) const { os << "unknown op"; return os; }
        };



        // first  -> original expression
        // second -> derivative expression &
        template <class T>
        using OriginalAndDerivativeExpression = std::pair<Expression<T>, DerivativeExpression<T> &>;
        template <class T>
        using OriginalAndDerivativeExpressionTable = std::vector<std::pair<Expression<T>, DerivativeExpression<T>>>;
        template <class T>
        inline OriginalAndDerivativeExpression<T> MakeOriginalAndDerivative(Expression<T> && input, DerivativeExpression<T> & dinput) {
            return std::pair<Expression<T>, DerivativeExpression<T> &>{input, dinput};
        }

        namespace {
            // functional op
            template <class OpTraitsT>
            class FunctionalOp : public OpBaseType<typename OpTraitsT::OutputType> {
            private:
                using InputIndices = typename OpTraitsT::InputIndices;

                using OutputType = typename OpTraitsT::OutputType;
                using InputTuple = typename OpTraitsT::InputTuple;
                static const int InputsNumber = OpTraitsT::InputsNumber;
                template <int I>
                struct InputTypeStruct {
                    using type = typename std::tuple_element<I, InputTuple>::type;
                };
                template <int I>
                using InputType = typename InputTypeStruct<I>::type;

            public:
                inline explicit FunctionalOp(OpTraitsT t) : OpBaseType<typename OpTraitsT::OutputType>(), _opTraits(t) {}

            private:
                template <int ...S>
                inline OutputType valueUsingSequence(std::vector<EHandle> && inputs, Sequence<S...>) const {
                    assert(inputs.size() == InputsNumber);
                    // inputs -> ResultType<inputs> -> output
                    return _opTraits.value(TraitsAboutOp<InputType<S>>::Result(graph->op(inputs[S])) ...);
                }

                template <int ...S>
                inline std::vector<EHandle> backPropagateGradientUsingSequence(
                    std::vector<EHandle> && inputs,
                    EHandle sumOfDOutputs,
                    Sequence<S...>) const {

                    // inputs + doutputs -> dinputs
                    std::tuple<DerivativeExpression<InputType<S>>...> dinputs;
                    _opTraits.derivatives(
                        graph->as<OutputType>(self), // output
                        graph->asDerived<OutputType>(sumOfDOutputs), // DOutput
                        MakeOriginalAndDerivative(graph->as<InputType<S>>(inputs[S]), std::get<S>(dinputs))...);
                    return std::vector<EHandle>{std::get<S>(dinputs).handle()...};
                }

            public:
                virtual std::ostream & toString(std::ostream & os) const {
                    return _opTraits.toString(os);
                }
                virtual OutputType value() const override {
                    return valueUsingSequence(graph->inputs(self), InputIndices());
                }
                virtual std::vector<EHandle> backPropagateGradient(EHandle sumOfDOutputs) const {
                    return backPropagateGradientUsingSequence(graph->inputs(self), sumOfDOutputs, InputIndices());
                }
                OpTraitsT _opTraits;
            };
        }

        // compose Expression using FunctionalOp
        template <class OpTraitsT, class InputT, class ...InputTs>
        inline Expression<typename OpTraitsT::OutputType> ComposeExpression(
            OpTraitsT opTraits,
            const Expression<InputT> & firstInput,
            const Expression<InputTs> & ... inputs){
            static_assert(OpTraitsT::InputsNumber == 1 + sizeof...(inputs), "inputs number mismatch!");
            return firstInput.g()->as<typename OpTraitsT::OutputType>(firstInput.g()->
                addNode(std::make_shared<FunctionalOp<OpTraitsT>>(opTraits), {
                firstInput.handle(), inputs.handle()...
            }));
        }

        // compose Expression (with no inputs) using FunctionalOp
        template <class OpTraitsT>
        inline Expression<typename OpTraitsT::OutputType> ComposeExpressionWithNoInputs(ExpressionGraph & graph,
            OpTraitsT opTraits){
            static_assert(OpTraitsT::InputsNumber == 0, "inputs number mismatch!");
            return graph.as<typename OpTraitsT::OutputType>(graph.addNode(std::make_shared<FunctionalOp<OpTraitsT>>(opTraits)));
        }








        // compose type/function
        namespace {
   
            template <class T, class ...ArgTs>
            struct TypeComposerTraits : public OpTraitsBase<T> {
                inline explicit TypeComposerTraits(const std::string & nm, ArgTs && ...args)
                    : args(std::forward_as_tuple(args...)), name(nm) {}
                inline T value() const {
                    return valueUsingSequence(typename SequenceGenerator<sizeof...(ArgTs)>::type());
                }
                template <int ...S>
                inline T && valueUsingSequence(Sequence<S...>) const {
                    return std::move(T(std::get<S>(args)...));
                }
                inline void derivatives(Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs) const {
                }
                virtual ostream & toString(ostream & os) const { os << name; return os; }
                std::tuple<ArgTs...> args;
                std::string name;
            };

            template <class FunctorT, class OutputT>
            struct FunctorComposerTraits : public OpTraitsBase<OutputT> {
                inline explicit FunctorComposerTraits(const std::string & nm, FunctorT && f) 
                    : fun(std::forward<FunctorT>(f)), name(nm) {}
                inline OutputType value() const {
                    return fun();
                }
                inline void derivatives(Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs) const {
                }
                virtual ostream & toString(ostream & os) const { os << name; return os; }
                FunctorT fun;
                std::string name;
            };

        }

        template <class T, class ...ArgTs>
        inline Expression<T> composeType(ExpressionGraph & graph, ArgTs && ... args) {
            return ComposeExpressionWithNoInputs(graph, TypeComposerTraits<T, ArgTs...>("aComposedType", std::forward<ArgTs>(args)...));
        }

        template <class T, class ...ArgTs>
        inline Expression<T> composeType(const std::string & name, ExpressionGraph & graph, ArgTs && ... args) {
            return ComposeExpressionWithNoInputs(graph, TypeComposerTraits<T, ArgTs...>(name, std::forward<ArgTs>(args)...));
        }

        template <class FunctorT>
        inline auto composeFunction(ExpressionGraph & graph, FunctorT && fun) -> Expression<decltype(fun())> {
            return ComposeExpressionWithNoInputs(graph, 
                FunctorComposerTraits<FunctorT, decltype(fun())>("aComposedFunction", std::forward<FunctorT>(fun)));
        }

        template <class FunctorT>
        inline Expression<std::result_of_t<FunctorT>> composeFunction(const std::string & name, ExpressionGraph & graph, FunctorT && fun) {
            return ComposeExpressionWithNoInputs(graph, 
                FunctorComposerTraits<FunctorT>(name, std::forward<FunctorT>(fun)));
        }


        // fill with scalar
        namespace  {
            template <class T>
            struct FillWithScalarTraits : public OpTraitsBase<DataStorageType<ResultType<T>>, T> {
                //static_assert(Fun_fillWithScalar_DefinedInDataTraits<T>::value, "fillWithScalar not defined!");
                using Scalar = DataScalarType<ResultType<T>>;
                inline explicit FillWithScalarTraits(const Scalar & ss) : s(ss) {}
                inline OutputType value(ResultType<T> t) const {
                    //return DataTraits<ResultType<T>>::fillWithScalar(std::move(t), s);
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
        inline Expression<DataStorageType<ResultType<T>>> fillWithScalar(Expression<T> t, DataScalarType<ResultType<T>> s) {
            return ComposeExpression(FillWithScalarTraits<T>(s), t);
        }



        // sum
        namespace  {

            template <class T>
            inline T SumAll(T && t){
                return t;
            }

            template <class T, class ... Ts>
            inline auto SumAll(T && t, Ts &&... ts) -> decltype(t + SumAll(ts...)) {
                return t + SumAll(ts...);
            }

            template <class ...T>
            struct SumResult;

            template <>
            struct SumResult<> {};

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
        inline Expression<SumResultType<Ts...>> generalSum(Expression<Ts>... inputs) {
            return ComposeExpression(SumTraits<Ts...>(), inputs...);
        }

        template <class T1, class T2>
        inline auto operator + (Expression<T1> a, Expression<T2> b) -> decltype(generalSum(a, b)) {
            return generalSum(a, b);
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



        // transpose
        namespace  {
            template <class T>
            struct TransposeResult {
                using type = DataStorageType<decltype(common::GeneralTranspose(std::declval<ResultType<T>>()))>;
            };
            template <class T>
            using TransposeResultType = typename TransposeResult<T>::type;

            template <class T>
            struct TransposeTraits : public OpTraitsBase<TransposeResultType<T>, T> {
                inline OutputType value(ResultType<T> from) const {
                    return common::GeneralTranspose(from);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> from) const {
                    from.second = transpose(sumOfDOutputs).eval().assign<DerivativeType<T>>();
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

            template <class T1, class T2>
            struct ProductResult{
                using type = DataStorageType<decltype(std::declval<ResultType<T1>>() * std::declval<ResultType<T2>>())>;
            };

            template <class T1, class T2>
            using ProductResultType = typename ProductResult<T1, T2>::type;

            // general product DataTraits
            template <class T1, class T2>
            struct ProductTraits : public OpTraitsBase<ProductResultType<T1, T2>, T1, T2> {
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
                virtual ostream & toString(ostream & os) const { os << "prod"; return os; }
            };

        }

        template <class T1, class T2>
        inline Expression<ProductResultType<T1, T2>> generalProd(const Expression<T1> & a, const Expression<T2> & b) {
            return ComposeExpression(ProductTraits<T1, T2>(), a, b);
        }

        template <class T1, class T2>
        inline auto operator * (const Expression<T1> & a, const Expression<T2> & b) -> decltype(generalProd(a, b)) {
            return generalProd(a, b);
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
                virtual ostream & toString(ostream & os) const { os << "subtract"; return os; }
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
                virtual ostream & toString(ostream & os) const { os << "multScalar[" << s << "]"; return os; }
                DataScalarType<ResultType<T>> s;
            };
        }

        template <class T>
        inline auto operator * (Expression<T> a, DataScalarType<ResultType<T>> b)
            ->Expression<decltype(std::declval<ResultType<T>>() * std::declval<DataScalarType<ResultType<T>>>())> {
            return ComposeExpression(MultScalar<T>(b), a);
        }

        template <class T>
        inline auto operator / (Expression<T> a, DataScalarType<ResultType<T>> b)
            ->Expression<decltype(std::declval<ResultType<T>>() * std::declval<DataScalarType<ResultType<T>>>())> {
            return ComposeExpression(MultScalar<T>(1.0/b), a);
        }

        template <class T>
        inline auto operator * (DataScalarType<ResultType<T>> b, Expression<T> a)
            ->Expression<decltype(std::declval<ResultType<T>>() * std::declval<DataScalarType<ResultType<T>>>())> {
            return ComposeExpression(MultScalar<T>(b), a);
        }

        // plus scalar
        namespace  {
            template <class T>
            struct PlusScalar : public OpTraitsBase<
                decltype(std::declval<ResultType<T>>() + std::declval<DataScalarType<ResultType<T>>>()), T> {

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
        inline auto operator + (Expression<T> a, DataScalarType<ResultType<T>> b)
            ->Expression<decltype(std::declval<ResultType<T>>() + std::declval<DataScalarType<ResultType<T>>>())> {
            return ComposeExpression(PlusScalar<T>(b), a);
        }

        template <class T>
        inline auto operator + (DataScalarType<ResultType<T>> b, Expression<T> a)
            ->Expression<decltype(std::declval<ResultType<T>>() + std::declval<DataScalarType<ResultType<T>>>())> {
            return ComposeExpression(PlusScalar<T>(b), a);
        }

        template <class T>
        inline auto operator - (Expression<T> a, DataScalarType<ResultType<T>> b)
            ->Expression<decltype(std::declval<ResultType<T>>() - std::declval<DataScalarType<ResultType<T>>>())> {
            return ComposeExpression(PlusScalar<T>(-b), a);
        }

        template <class T>
        inline auto operator - (DataScalarType<ResultType<T>> b, Expression<T> a) -> decltype((-a) + b) {
            return (-a) + b;
        }



        // cwise product
        namespace {

            template <class T1, class T2>
            using CWiseProductResultType = decltype(common::CWiseProd(std::declval<T1>(), std::declval<T2>()));


            // cwise product traits
            template <class T1, class T2>
            struct CWiseProductTraits : public OpTraitsBase<CWiseProductResultType<T1, T2>, T1, T2> {
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
        inline Expression<CWiseProductResultType<T1, T2>> cwiseProd(const Expression<T1> & a, const Expression<T2> & b) {
            return ComposeExpression(CWiseProductTraits<T1, T2>(), a, b);
        }



        // pow
        namespace  {

            using std::pow;

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
                    input.second = cwiseProd(pow(input.first, exponent - 1.0) * exponent, sumOfDOutputs).eval();
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
        inline auto CWiseQuotient(const Expression<T1> & a, const Expression<T2> & b) -> decltype(cwiseProd(a, pow(b, -1))) {
            return cwiseProd(a, pow(b, -1));
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
                    input.second = cwiseProd(sumOfDOutputs, output).eval();
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


        // sigmoid
        namespace {
            
            using std::exp;
            
            template <class T>
            struct SigmoidResult {
                using ScalarType = DataScalarType<T>;
                using type = decltype(std::declval<ScalarType>() / 
                    (std::declval<ScalarType>() + exp(-std::declval<ResultType<T>>())));
            };
            template <class T>
            using SigmoidResultType = typename SigmoidResult<T>::type;

            template <class T>
            struct SigmoidTraits : public OpTraitsBase<SigmoidResultType<T>, T> {
                inline OutputType value(ResultType<T> input) const {
                    return 1.0 / (1.0 + exp(-input));
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
        inline Expression<SigmoidResultType<T>> sigmoid(const Expression<T> & e) {
            return ComposeExpression(SigmoidTraits<T>(), e);
        }

        
        // cwiseSelect
        namespace {

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectResult {
                using type = decltype(common::CWiseSelect(
                    std::declval<ResultType<CondT>>(), 
                    std::declval<ResultType<IfT>>(), 
                    std::declval<ResultType<ElseT>>()));
            };
            template <class CondT, class IfT, class ElseT>
            using CWiseSelectResultType = typename CWiseSelectResult<CondT, IfT, ElseT>::type;

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectTraits : public OpTraitsBase<CWiseSelectResultType<CondT, IfT, ElseT>, CondT, IfT, ElseT> {
                inline OutputType value(ResultType<CondT> cond, ResultType<IfT> ifval, ResultType<ElseT> elseval) const {
                    return common::CWiseSelect(cond, ifval, elseval);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<CondT> cond,
                    OriginalAndDerivativeExpression<IfT> ifval,
                    OriginalAndDerivativeExpression<ElseT> elseval) const {
                    ifval.second = cwiseSelect(cond.first, sumOfDOutputs, 0).eval();
                    elseval.second = cwiseSelect(cond.first, 0, sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "cwiseSelect"; return os; }
            };

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectIfRetConstTraits : public OpTraitsBase<CWiseSelectResultType<CondT, IfT, ElseT>, CondT, ElseT> {
                inline explicit CWiseSelectIfRetConstTraits(IfT ifv) : ifval(ifv) {}
                inline OutputType value(ResultType<CondT> cond, ResultType<ElseT> elseval) const {
                    return common::CWiseSelect(cond, ifval, elseval);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<CondT> cond,
                    OriginalAndDerivativeExpression<ElseT> elseval) const {
                    elseval.second = cwiseSelect(cond.first, 0, sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "cwiseSelect[ifthen:" << ifval << "]"; return os; }
                IfT ifval;
            };

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectElseRetConstTraits : public OpTraitsBase<CWiseSelectResultType<CondT, IfT, ElseT>, CondT, IfT> {
                inline explicit CWiseSelectElseRetConstTraits(ElseT elsev) : elseval(elsev) {}
                inline OutputType value(ResultType<CondT> cond, ResultType<IfT> ifval) const {
                    return common::CWiseSelect(cond, ifval, elseval);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<CondT> cond,
                    OriginalAndDerivativeExpression<IfT> ifval) const {
                    ifval.second = cwiseSelect(cond.first, sumOfDOutputs, 0).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "cwiseSelect[elsethen:" << elseval << "]"; return os; }
                ElseT elseval;
            };
        }

        template <class CondT, class IfT, class ElseT>
        inline Expression<CWiseSelectResultType<CondT, IfT, ElseT>> cwiseSelect(
            const Expression<CondT> & cond, const Expression<IfT> & ifval, const Expression<ElseT> & elseval) {
            return ComposeExpression(CWiseSelectTraits<CondT, IfT, ElseT>(), cond, ifval, elseval);
        }

        template <class CondT, class ElseT>
        inline Expression<CWiseSelectResultType<CondT, DataScalarType<ElseT>, ElseT>> cwiseSelect(
            const Expression<CondT> & cond, DataScalarType<ElseT> ifval, const Expression<ElseT> & elseval) {
            return ComposeExpression(CWiseSelectIfRetConstTraits<CondT, DataScalarType<ElseT>, ElseT>(
                static_cast<DataScalarType<ElseT>>(ifval)), cond, elseval);
        }

        template <class CondT, class IfT>
        inline Expression<CWiseSelectResultType<CondT, IfT, DataScalarType<IfT>>> cwiseSelect(
            const Expression<CondT> & cond, const Expression<IfT> & ifval, DataScalarType<IfT> elseval) {
            return ComposeExpression(CWiseSelectElseRetConstTraits<CondT, IfT, DataScalarType<IfT>>(
                static_cast<DataScalarType<IfT>>(elseval)), cond, ifval);
        }


        // abs
        namespace {

            using std::abs;

            template <class T>
            struct AbsResult {
                using type = decltype(abs(std::declval<T>()));
            };
            template <class T>
            using AbsResultType = typename AbsResult<T>::type;

            template <class T>
            struct AbsTraits : public OpTraitsBase<AbsResultType<T>, T> {
                inline OutputType value(ResultType<T> input) const {
                    return abs(input);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input) const {
                    input.second = cwiseSelect(input.first, sumOfDOutputs, -sumOfDOutputs);
                }
                virtual ostream & toString(ostream & os) const { os << "abs"; return os; }
            };
        }

        template <class T>
        inline Expression<AbsResultType<T>> abs(const Expression<T> & e) {
            return ComposeExpression(AbsTraits<T>(), e);
        }



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


        // for eigen
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
                virtual ostream & toString(ostream & os) const { os << "matrixToArray"; return os; }
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
                virtual ostream & toString(ostream & os) const { os << "arrayToMatrix"; return os; }
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


        //// conjugate twin operators
        //namespace {
        //    template <class To, class From, class FunForwardT, class FunBackwardT>
        //    struct TwinTraits :
        //        public OpTraitsBase<To, From> {
        //        inline TwinTraits(FunForwardT funa, FunBackwardT funb) : forwardFun(funa), backwardFun(funb) {}
        //        inline OutputType value(ResultType<From> in) const {
        //            return forwardFun(in);
        //        }
        //        inline void derivatives(
        //            Expression<OutputType> output,
        //            DerivativeExpression<OutputType> sumOfDOutputs,
        //            OriginalAndDerivativeExpression<InputTypes<0>> input) const {
        //            input.second =
        //                ComposeExpression(TwinTraits<From, To, FunBackwardT, FunForwardT>(backwardFun, forwardFun), 
        //                sumOfDOutputs);
        //        }
        //        FunForwardT forwardFun;
        //        FunBackwardT backwardFun;
        //    };
        //}


        //template <class To, class From, class FunForwardT, class FunBackwardT>
        //inline Expression<To> composeTwinOp(const Expression<From> & e, FunForwardT forwardFun, FunBackwardT backwardFun){
        //    return ComposeExpression(TwinTraits<To, From, FunForwardT, FunBackwardT>(forwardFun, backwardFun), e);
        //}

       
        

    }
}
 
#endif