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
                //static_assert(Fun_fillWithScalar_DefinedInDataTraits<T>::value, "fillWithScalar not defined!");
                using Scalar = DataScalarType<T>;
                inline explicit FillWithScalarWhenScalarIsConstTraits(const Scalar & ss) : s(ss) {}
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
        inline Expression<DataStorageType<T>> fillWithScalar(Expression<T> t, Expression<DataScalarType<T>> & s) {
            return ComposeExpression(FillWithScalarTraits<T>(), t, s);
        }

        template <class T>
        inline Expression<DataStorageType<T>> fillWithScalar(Expression<T> t, DataScalarType<T> s) {
            return ComposeExpression(FillWithScalarWhenScalarIsConstTraits<T>(s), t);
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
        inline auto generalSum(Expression<Ts>... inputs) 
            -> decltype(ComposeExpression(SumTraits<Ts...>(), inputs...)) {
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
            struct TransposeTraits 
                : public OpTraitsBase<decltype(common::Eval(common::GeneralTranspose(std::declval<ResultType<T>>()))), T> {
                inline OutputType value(ResultType<T> from) const {
                    return common::Eval(common::GeneralTranspose(from));
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
        inline auto transpose(const Expression<T> & e) -> decltype(ComposeExpression(TransposeTraits<T>(), e)) {
            return ComposeExpression(TransposeTraits<T>(), e);
        }




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
        inline auto cwiseProd(const Expression<T1> & a, const Expression<T2> & b) -> decltype(ComposeExpression(CWiseProductTraits<T1, T2>(), a, b)) {
            return ComposeExpression(CWiseProductTraits<T1, T2>(), a, b);
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
                    input.second = cwiseProd(cwiseProd(1.0 - output, 1.0 - output), sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "tanh"; return os; }
            };
        }

        template <class T>
        inline auto tanh(const Expression<T> & e) -> decltype(ComposeExpression(TanhTraits<T>(), e)) {
            return ComposeExpression(TanhTraits<T>(), e);
        }


        
        // cwiseSelect
        namespace {

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectTraitsWithoutConsts 
                : public OpTraitsBase<decltype(common::CWiseSelect(
                std::declval<ResultType<CondT>>(),
                std::declval<ResultType<IfT>>(),
                std::declval<ResultType<ElseT>>())), CondT, IfT, ElseT> {

                inline OutputType value(ResultType<CondT> cond, ResultType<IfT> ifval, ResultType<ElseT> elseval) const {
                    return common::CWiseSelect(cond, ifval, elseval);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<CondT> cond,
                    OriginalAndDerivativeExpression<IfT> ifval,
                    OriginalAndDerivativeExpression<ElseT> elseval) const {
                    ifval.second = cwiseSelect(cond.first, sumOfDOutputs, 0.0).eval();
                    elseval.second = cwiseSelect(cond.first, 0.0, sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "cwiseSelect"; return os; }
            };

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectWhenIfRetIsConstTraits 
                : public OpTraitsBase<decltype(common::CWiseSelect(
                std::declval<ResultType<CondT>>(),
                std::declval<ResultType<DataStorageType<IfT>>>(),
                std::declval<ResultType<ElseT>>())), CondT, ElseT> {

                inline explicit CWiseSelectWhenIfRetIsConstTraits(const IfT & ifv) : ifval(ifv) {}
                inline OutputType value(ResultType<CondT> cond, ResultType<ElseT> elseval) const {
                    return common::CWiseSelect(cond, ifval, elseval);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<CondT> cond,
                    OriginalAndDerivativeExpression<ElseT> elseval) const {
                    elseval.second = cwiseSelect(cond.first, 0.0, sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "cwiseSelect[ifthen:" << ifval << "]"; return os; }
                DataStorageType<IfT> ifval;
            };

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectWhenElseRetIsConstTraits 
                : public OpTraitsBase<decltype(common::CWiseSelect(
                std::declval<ResultType<CondT>>(),
                std::declval<ResultType<IfT>>(),
                std::declval<ResultType<DataStorageType<ElseT>>>())), CondT, IfT> {

                inline explicit CWiseSelectWhenElseRetIsConstTraits(const ElseT & elsev) : elseval(elsev) {}
                inline OutputType value(ResultType<CondT> cond, ResultType<IfT> ifval) const {
                    return common::CWiseSelect(cond, ifval, elseval);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<CondT> cond,
                    OriginalAndDerivativeExpression<IfT> ifval) const {
                    ifval.second = cwiseSelect(cond.first, sumOfDOutputs, 0.0).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "cwiseSelect[elsethen:" << elseval << "]"; return os; }
                DataStorageType<ElseT> elseval;
            };

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectTraits {
                inline static Expression<IfT> Compose(const CondT & cond, const IfT & ifval, const ElseT & elseval) {
                    SHOULD_NEVER_BE_INSTANCIATED();
                }
            };

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectTraits<Expression<CondT>, IfT, Expression<ElseT>> {
                inline static auto Compose(const Expression<CondT> & cond, const IfT & ifval, const Expression<ElseT> & elseval)
                -> decltype(ComposeExpression(CWiseSelectWhenIfRetIsConstTraits<CondT, IfT, ElseT>(ifval), cond, elseval)) {
                    return ComposeExpression(CWiseSelectWhenIfRetIsConstTraits<CondT, IfT, ElseT>(ifval), cond, elseval);
                }
            };

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectTraits<Expression<CondT>, Expression<IfT>, ElseT> {
                inline static auto Compose(const Expression<CondT> & cond, const Expression<IfT> & ifval, const ElseT & elseval)
                -> decltype(ComposeExpression(CWiseSelectWhenElseRetIsConstTraits<CondT, IfT, ElseT>(elseval), cond, ifval)) {
                    return ComposeExpression(CWiseSelectWhenElseRetIsConstTraits<CondT, IfT, ElseT>(elseval), cond, ifval);
                }
            };

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectTraits<Expression<CondT>, Expression<IfT>, Expression<ElseT>> {
                inline static auto Compose(const Expression<CondT> & cond, const Expression<IfT> & ifval, const Expression<ElseT> & elseval)
                -> decltype(ComposeExpression(CWiseSelectTraitsWithoutConsts<CondT, IfT, ElseT>(), cond, ifval, elseval)) {
                    return ComposeExpression(CWiseSelectTraitsWithoutConsts<CondT, IfT, ElseT>(), cond, ifval, elseval);
                }
            };

        }

        template <class CondT, class IfT, class ElseT>
        inline auto cwiseSelect(const CondT & cond, const IfT & ifval, const ElseT & elseval)
            -> decltype(CWiseSelectTraits<CondT, IfT, ElseT>::Compose(cond, ifval, elseval)) {
            return CWiseSelectTraits<CondT, IfT, ElseT>::Compose(cond, ifval, elseval);
        }

        /*template <class CondT, class IfT, class ElseT>
        inline auto cwiseSelect(const Expression<CondT> & cond, const Expression<IfT> & ifval, const Expression<ElseT> & elseval) 
            -> decltype(ComposeExpression(CWiseSelectTraitsWithoutConsts<CondT, IfT, ElseT>(), cond, ifval, elseval)) {
            return ComposeExpression(CWiseSelectTraitsWithoutConsts<CondT, IfT, ElseT>(), cond, ifval, elseval);
        }

        template <class CondT, class ElseT, class IfT>
        inline auto cwiseSelect(const Expression<CondT> & cond, const IfT & ifval, const Expression<ElseT> & elseval, 
            std::enable_if_t<!IsExpression<IfT>::value, int> = 0)
            -> decltype(ComposeExpression(CWiseSelectWhenIfRetIsConstTraits<CondT, IfT, ElseT>(ifval), cond, elseval)) {
            return ComposeExpression(CWiseSelectWhenIfRetIsConstTraits<CondT, IfT, ElseT>(ifval), cond, elseval);
        }

        template <class CondT, class IfT, class ElseT>
        inline auto cwiseSelect(const Expression<CondT> & cond, const Expression<IfT> & ifval, const ElseT & elseval, 
            std::enable_if_t<!IsExpression<IfT>::value, int> = 0)
            -> decltype(ComposeExpression(CWiseSelectWhenElseRetIsConstTraits<CondT, IfT, ElseT>(elseval), cond, ifval)) {
            return ComposeExpression(CWiseSelectWhenElseRetIsConstTraits<CondT, IfT, ElseT>(elseval), cond, ifval);
        }*/




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

        // mapping
        namespace {
            template <class To, class From, class MappingFunT, class ReverseMappingFunT>
            struct MappingTraits : public OpTraitsBase<To, From> {
                static_assert(IsStorageType<To>::value && IsStorageType<From>::value, "To and From must both be storage types!");
                inline explicit MappingTraits(MappingFunT mf, ReverseMappingFunT rmf, 
                    const std::string & nm = "map", const std::string & rnm = "map") 
                    : mappingFun(mf), reverseMappingFun(rmf), name(nm), reversedName(rnm) {}
                inline To value(ResultType<From> from) const {
                    return mappingFun(from);
                }
                inline void derivatives(
                    Expression<To> to,
                    DerivativeExpression<To> sumOfDOutputs,
                    OriginalAndDerivativeExpression<From> from) const {
                    from.second = ComposeExpression(MappingTraits<From, To, ReverseMappingFunT, MappingFunT>
                        (reverseMappingFun, mappingFun), sumOfDOutputs);
                }
                virtual ostream & toString(ostream & os) const { os << name; return os; }
                MappingFunT mappingFun;
                ReverseMappingFunT reverseMappingFun;
                std::string name, reversedName;
            };
        }

        template <class To, class From, class MappingFunT, class ReverseMappingFunT>
        inline Expression<To> composeMappingFunction(const Expression<From> & from, 
            MappingFunT mf, ReverseMappingFunT rmf, const std::string & nm = "map", const std::string & rnm = "map") {
            return ComposeExpression(MappingTraits<To, From, MappingFunT, ReverseMappingFunT>(mf, rmf, nm, rnm), from);
        }
        
        

    }
}
 
#endif