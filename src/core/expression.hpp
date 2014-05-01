#ifndef PANORAMIX_CORE_EXPRESSION_HPP
#define PANORAMIX_CORE_EXPRESSION_HPP

#include <iostream>
#include <utility>
#include <algorithm>
#include <numeric>
#include <initializer_list>

#include "template_utilities.hpp"
#include "mesh.hpp"
#include "data_traits.hpp"

namespace panoramix {
    namespace core {

        class ExpressionGraph;
        struct Op;

        using GraphType = Mesh<std::shared_ptr<Op>>;
        using EHandle = GraphType::VertHandle;
        using CHandle = GraphType::HalfHandle;
        using std::ostream;


        template <class T> struct OpWithValue;
        template <class T> struct OpWithCache;
        template <class T> struct Expression;
        
        // operator base
        struct Op {

            // the graph
            ExpressionGraph * graph;

            // the handle of self node
            EHandle self;

            // stream to string
            virtual std::ostream & toString(std::ostream & os) const = 0;


            // execute the operator
            // called in forward propagation
            virtual void forwardPropagateExecute() {}


            // make input derivative expressions based on inputs and output derivative expressions
            // the input derivatives should be calculated based on 
            //  1. the sum of outputDerivs 
            //  2. the inputs
            // returns
            //  1. the sum of outputDerivs
            //  2. derivative Ops corresponding to inputs
            // dinputs size MUST be same with inputs number
            virtual EHandle homomorphicSum(const std::vector<EHandle> & doutputs) const = 0;
            virtual std::vector<EHandle> backPropagateGradient(EHandle sumOfDOutputs) const = 0;


            // make a new Op which makes a scalar one
            // requires that the Op represents a scalar
            virtual std::shared_ptr<Op> one() const = 0;


            // check whether this op has the value type
            template <class T>
            bool hasValue() const { 
                return dynamic_cast<const OpWithValue<T> *>(this) != nullptr; 
            }

            // convert to op with a value to retrieve a value
            template <class T>
            T value() const {
                return static_cast<const OpWithValue<T> *>(this)->value();
            }

            // check whether has cache
            template <class T>
            bool hasCache() const {
                return dynamic_cast<const OpWithCache<T> *>(this) != nullptr;
            }

            // convert to op with a cache to retrieve cache
            template <class T>
            const storage_type<T> & cache() const {
                return static_cast<const OpWithCache<T> *>(this)->cache;
            }

        };

        // base op for expressions
        template <class T>
        struct OpWithValue : public Op {

            // compute value
            virtual T value() const = 0;

            // do nothing
            virtual void forwardPropagateExecute() {}

            // make a scalar one to initialize automatic derivation
            virtual std::shared_ptr<Op> one() const override {
                static std::shared_ptr<Op> _one = 
                    std::make_shared<OpWithConstant<DerivativeType<T>>>(traits<DerivativeType<T>>::one());
                return _one; 
            }

            // homomorphicSum
            virtual EHandle homomorphicSum(const std::vector<EHandle> & doutputs) const override {
                assert(doutputs.size() > 0 && "doutputs size is zero!");
                if (doutputs.size() == 1)
                    return doutputs.front();
                return HSum<DerivativeType<T>>(graph, doutputs);
            }
        };        

        // using a value cache to store data
        template <class T>
        struct OpWithCache : public OpWithValue<T> {
            using ST = storage_type<T>;
            static_assert(std::is_assignable<ST&, T>::value, "T cannot be assigned to ST&!");
            inline OpWithCache(const ST & c = ST(), const std::string nm = "")
                : cache(c), givenName(nm)
            {}

            virtual std::ostream & toString(std::ostream & os) const override { 
                os << givenName; 
                return os; 
            }

            // compute value and store in cache
            virtual void forwardPropagateExecute() { cache = traits<T>::eval(value()); }

            ST cache;
            std::string givenName;
        };

        // constant
        template <class T>
        struct OpWithConstant : public OpWithCache<T> {
            static_assert(is_storage_type<T>::value, "T must be a storage type!");
            inline OpWithConstant(const T & c, const std::string nm = "") : OpWithCache<T>(c, nm) {}
            virtual std::ostream & toString(std::ostream & os) const override {
                os << givenName << "[" << cache << "]";
                return os;
            }
            virtual T value() const override { return cache; }
            virtual void forwardPropagateExecute() override {}
            virtual std::vector<EHandle> backPropagateGradient(EHandle sumOfDOutputs) const override { return std::vector<EHandle>(); }
        };

        // using a reference type to reference data
        template <class T>
        struct OpWithReference : public OpWithValue<const T &> {
            inline OpWithReference(const T & c, const std::string & nm = "ref") 
                : OpWithValue<const T &>(), ref(c), givenName(nm) {}
            
            virtual ostream & toString(ostream & os) const override {
                os << givenName;
                return os;
            }

            virtual const T & value() const override { return ref; }
            virtual void forwardPropagateExecute() override {}
            virtual std::vector<EHandle> backPropagateGradient(EHandle sumOfDOutputs) const override { return std::vector<EHandle>(); }

            const T & ref;
            std::string givenName;
        };










        // traits about op
        struct result_retrieved_by_value {};
        struct result_retrieved_by_cache {};
        template <class T>
        struct result_tag {
            using type = 
            std::conditional_t < std::is_lvalue_reference<T>::value, result_retrieved_by_value,
            std::conditional_t < traits<T>::should_be_cached, result_retrieved_by_cache,
            result_retrieved_by_value >> ;
        };

        template <class T, class ResultTag = typename result_tag<T>::type>
        struct TraitsAboutOp {};

        template <class T> // by value
        struct TraitsAboutOp<T, result_retrieved_by_value> { 
            using OpType = OpWithValue<T>;
            using ResultType = T;
            static T Result(Op const & op) { 
                if (!op.hasValue<T>()){
                    throw std::runtime_error(std::string("it is determined that ") + typeid(T).name() + 
                        " is retrieved by value(), but current op doesn't has the correct value!");
                }
                return op.value<T>(); 
            }
        };

        template <class T>
        struct TraitsAboutOp<T, result_retrieved_by_cache> { // by cache
            using OpType = OpWithCache<T>;
            using ResultType = const storage_type<T> &;
            static ResultType Result(Op const & op) { 
                if (!op.hasCache<T>()){
                    throw std::runtime_error(std::string("it is determined that ") + typeid(T).name() +
                        " is retrieved by cache(), but current op doesn't has the correct cache!");
                }
                return op.cache<T>(); 
            }
        };


        // smartly determines whether cache should be used
        template <class T>
        using OpBaseType = typename TraitsAboutOp<T>::OpType;
        
        // smartly determins the result return type
        template <class T>
        using ResultType = typename TraitsAboutOp<T>::ResultType;

        // type of derivative
        template <class T>
        using DerivativeType = storage_type<T>;

        // type of derivative expression
        template <class T>
        using DerivativeExpression = Expression<storage_type<T>>;
        
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
        // second -> derivative expression 
        template <class T>
        using OriginalAndDerivativeExpression = std::pair<Expression<T>, DerivativeExpression<T> &>;
        template <class T>
        OriginalAndDerivativeExpression<T>
            MakeOriginalAndDerivative(Expression<T> && input, DerivativeExpression<T> & dinput) {
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
                        graph->asDerivative<OutputType>(sumOfDOutputs), // DOutput
                        MakeOriginalAndDerivative(graph->as<InputType<S>>(inputs[S]), std::get<S>(dinputs))...);
                    return std::vector<EHandle>{std::get<S>(dinputs).handle()...};
                }

            public:
                // get name
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
















        // op traits
        // expression cast
        namespace {
            // implicit assign // =
            template <class A, class B>
            struct ExpressionAssignTraits : public OpTraitsBase<A, B> {
                inline A value(ResultType<B> b) const {
                    return b;
                }
                inline void derivatives(
                    Expression<A> output,
                    DerivativeExpression<A> sumOfDOutputs,
                    OriginalAndDerivativeExpression<B> from) const {
                    from.second =
                        ComposeExpression(ExpressionAssignTraits<DerivativeType<B>, DerivativeType<A>>(), sumOfDOutputs);
                }
                virtual ostream & toString(ostream & os) const { os << "assign"; return os; }
            };

            // cast // data_cast
            template <class To, class From>
            struct ExpressionCastTraits : public OpTraitsBase<To, From> {
                inline To value(ResultType<From> from) const {
                    return data_cast<To>(from);
                }
                inline void derivatives(
                    Expression<To> output,
                    DerivativeExpression<To> sumOfDOutputs,
                    OriginalAndDerivativeExpression<From> from) const {
                    from.second = 
                        ComposeExpression(ExpressionCastTraits<DerivativeType<From>, DerivativeType<To>>(), 
                        sumOfDOutputs);
                }
                virtual ostream & toString(ostream & os) const { os << "cast"; return os; }
            };

            // eval // .eval()
            template <class T>
            struct ExpressionEvalTraits : public OpTraitsBase<storage_type<ResultType<T>>, T> {
                inline OutputType value(ResultType<T> from) const {
                    return traits<ResultType<T>>::eval(from);
                }
                static_assert(std::is_same<DerivativeExpression<OutputType>,
                    DerivativeExpression<T>>::value, "should be the same!");
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> from) const {
                    from.second = sumOfDOutputs;
                }
                virtual ostream & toString(ostream & os) const { os << "eval"; return os; }
            };
        }

        template <class To, class From>
        inline std::enable_if_t<std::is_same<To, From>::value, Expression<To>> 
            expression_assign_to(const Expression<From> & from) {
            return from;
        }

        template <class To, class From>
        inline std::enable_if_t<!std::is_same<To, From>::value, Expression<To>> 
            expression_assign_to(const Expression<From> & from) {
            return ComposeExpression(ExpressionAssignTraits<To, From>(), from);
        }

        template <class To, class From>
        inline std::enable_if_t<std::is_same<To, From>::value, Expression<To>> 
            expression_cast(const Expression<From> & from) {
            return from;
        }

        template <class To, class From>
        inline std::enable_if_t<!std::is_same<To, From>::value, Expression<To>> 
            expression_cast(const Expression<From> & from) {
            return ComposeExpression(ExpressionCastTraits<To, From>(), from);
        }

        template <class T>
        inline std::enable_if_t<is_storage_type<ResultType<T>>::value,
            typename ExpressionEvalTraits<T>::OutputExpressionType> 
            expression_eval(const Expression<T> & from) {
            return from;
        }

        template <class T>
        inline std::enable_if_t<!is_storage_type<ResultType<T>>::value, 
            typename ExpressionEvalTraits<T>::OutputExpressionType>
            expression_eval(const Expression<T> & from) {
            return ComposeExpression(ExpressionEvalTraits<T>(), from);
        }


        // set constant
        namespace {
            template <class T>
            struct SetConstantTraits : public OpTraitsBase<storage_type<ResultType<T>>, T> {
                using Scalar = scalar_type<ResultType<T>>;
                inline explicit SetConstantTraits(const Scalar & ss) : s(ss) {}
                inline OutputType value(ResultType<T> t) const {
                    return traits<ResultType<T>>::fill(t, s);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> from) const {
                    from.second = DerivativeExpression<T>(); // disconnected
                }
                virtual ostream & toString(ostream & os) const { os << "setConstant[" << s << "]"; return os; }
                Scalar s;
            };
        }

        template <class T>
        inline typename SetConstantTraits<T>::OutputExpressionType set_constant(Expression<T> t, scalar_type<ResultType<T>> s) {
            return ComposeExpression(SetConstantTraits<T>(s), t);
        }


        // linear combination
        namespace {

            template <class T>
            inline T SumAll(T && t){
                return t;
            }

            template <class T, class ... Ts>
            inline auto SumAll(T && t, Ts &&... ts) 
                -> decltype(t + SumAll(ts...)) {
                return t + SumAll(ts...);
            }

            template <class ...Ts>
            struct LinearCombResult {
                using type = decltype(SumAll(std::declval<ResultType<Ts>>() * std::declval<scalar_type<ResultType<Ts>>>() ...));
            };
            template <class ...Ts>
            using LinearCombResultType = typename LinearCombResult<Ts...>::type;

            template <class ...Ts>
            struct LinearCombTraits : public OpTraitsBase<LinearCombResultType<Ts...>, Ts...> {
                inline explicit LinearCombTraits(scalar_type<ResultType<Ts>> ... cs) : coeffs(std::make_tuple(cs...)) {}
                inline OutputType value(ResultType<Ts> ... inputs) const {
                    return valueUsingSequence(typename SequenceGenerator<InputsNumber>::type(), inputs...);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<Ts> ... inputs) const {
                    derivativesUsingSequence(typename SequenceGenerator<InputsNumber>::type(), sumOfDOutputs, inputs...);
                }

                virtual ostream & toString(ostream & os) const { os << "linsum"; return os; }

            private:
                template <int ...S>
                inline OutputType valueUsingSequence(Sequence<S...>, ResultType<Ts> ... inputs) const {
                    return SumAll((inputs * std::get<S>(coeffs)) ...);
                }

                template <int ...S>
                inline void derivativesUsingSequence(Sequence<S...>,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<Ts> ... inputs) const {
                    std::tie(inputs.second ...) = 
                        std::tie((sumOfDOutputs * static_cast<scalar_type<DerivativeType<OutputType>>>(std::get<S>(coeffs)))
                        .castTo<DerivativeType<Ts>>()
                        .eval()...);
                }

            public:
                std::tuple<scalar_type<ResultType<Ts>>...> coeffs;
            };

           

        }

        template <class T1, class T2>
        inline typename LinearCombTraits<T1, T2>::OutputExpressionType operator + (Expression<T1> a, Expression<T2> b) {
            return ComposeExpression(LinearCombTraits<T1, T2>(1.0, 1.0), a, b);
        }
        
        template <class T1, class T2>
        inline typename LinearCombTraits<T1, T2>::OutputExpressionType operator - (Expression<T1> a, Expression<T2> b) {
            return ComposeExpression(LinearCombTraits<T1, T2>(1.0, -1.0), a, b);
        }

        template <class T>
        inline typename LinearCombTraits<T>::OutputExpressionType operator - (Expression<T> a) {
            return ComposeExpression(LinearCombTraits<T>(-1.0), a);
        }

        template <class T>
        inline typename LinearCombTraits<T>::OutputExpressionType operator * (Expression<T> a, scalar_type<T> b) {
            return ComposeExpression(LinearCombTraits<T>(b), a);
        }

        template <class T>
        inline typename LinearCombTraits<T>::OutputExpressionType operator / (Expression<T> a, scalar_type<T> b) {
            return ComposeExpression(LinearCombTraits<T>(1.0/b), a);
        }

        template <class T>
        inline typename LinearCombTraits<T>::OutputExpressionType operator * (scalar_type<T> b, Expression<T> a) {
            return ComposeExpression(LinearCombTraits<T>(b), a);
        }

        namespace {

            template <class T, class EHandleIteratorT>
            inline EHandle HSum2(ExpressionGraph * graph, EHandleIteratorT inpuths) {
                return ComposeExpression(LinearCombTraits<T, T>(1, 1),
                    graph->as<T>(*(inpuths)),
                    graph->as<T>(*(inpuths+1))).eval().handle();
            }
            template <class T, class EHandleIteratorT>
            inline EHandle HSum3(ExpressionGraph * graph, EHandleIteratorT inpuths) {
                return ComposeExpression(LinearCombTraits<T, T, T>(1, 1, 1),
                    graph->as<T>(*(inpuths)),
                    graph->as<T>(*(inpuths+1)),
                    graph->as<T>(*(inpuths+2))).eval().handle();
            }
            template <class T, class EHandleIteratorT>
            inline EHandle HSum4(ExpressionGraph * graph, EHandleIteratorT inpuths) {
                return ComposeExpression(LinearCombTraits<T, T, T, T>(1, 1, 1, 1),
                    graph->as<T>(*(inpuths)),
                    graph->as<T>(*(inpuths+1)),
                    graph->as<T>(*(inpuths+2)),
                    graph->as<T>(*(inpuths+3))).eval().handle();
            }
            template <class T, class EHandleIteratorT>
            inline EHandle HSum5(ExpressionGraph * graph, EHandleIteratorT inpuths) {
                return ComposeExpression(LinearCombTraits<T, T, T, T, T>(1, 1, 1, 1, 1),
                    graph->as<T>(*(inpuths)),
                    graph->as<T>(*(inpuths+1)),
                    graph->as<T>(*(inpuths+2)),
                    graph->as<T>(*(inpuths+3)),
                    graph->as<T>(*(inpuths+4))).eval().handle();
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


        // sum all elements
        namespace {
            template <class T>
            struct SumElementsTraits : public OpTraitsBase<scalar_type<ResultType<T>>, T> {
                inline OutputType value(ResultType<T> t) const {
                    return t.sum(); // only eigen is supported
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> from) const {
                    from.second = set_constant(from.first, 1);
                }
                virtual ostream & toString(ostream & os) const { os << "sumElements"; return os; }
            };
        }

        template <class T>
        Expression<scalar_type<ResultType<T>>> sum_elements(const Expression<T> & e){
            return ComposeExpression(SumElementsTraits<T>(), e);
        }


        // transpose
        namespace {
            template <class T>
            struct TransposeTraits : public OpTraitsBase<decltype(traits<ResultType<T>>::transpose(std::declval<ResultType<T>>())), T> {
                inline OutputType value(ResultType<T> from) const {
                    return traits<ResultType<T>>::transpose(from);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> from) const {
                    from.second = transpose(sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "transpose"; return os; }
            };

            template <class T>
            using TransposeResultType = typename TransposeTraits<T>::OutputType;
        }

        template <class T>
        Expression<TransposeResultType<T>> 
            transpose(const Expression<T> & e) {
            return ComposeExpression(TransposeTraits<T>(), e);
        }


        // product
        namespace {

            // cwise product
            template <class T>
            inline T CWiseProdAll(T && t){
                return t;
            }

            template <class T, class ... Ts>
            inline auto CWiseProdAll(T && t, Ts &&... ts)
                -> decltype(data_cwise_product(t, CWiseProdAll(ts...))) {
                return data_cwise_product(t, CWiseProdAll(ts...));
            }


            template <class T, class ...Ts>
            struct CWiseProductionResult {
                using ScalarType = scalar_type<ResultType<T>>;
                static_assert(Sequence<std::is_same<scalar_type<ResultType<Ts>>, ScalarType>::value...>::All,
                     "all scalar types of Ts MUST be the same!");
                using type = 
                    decltype(CWiseProdAll(std::declval<ResultType<T>>(), std::declval<ResultType<Ts>>() ...));
            };
            template <class T, class ...Ts>
            using CWiseProductionResultType = typename CWiseProductionResult<T, Ts...>::type;
                

            // cwise product traits
            template <class T, class ...Ts>
            struct CWiseProductionTraits : public OpTraitsBase<CWiseProductionResultType<T, Ts...>, T, Ts...> {
                inline OutputType value(ResultType<T> input1, ResultType<Ts> ... inputs) const {
                    return CWiseProdAll(input1, inputs...);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input1,
                    OriginalAndDerivativeExpression<Ts> ... inputs) const {
                    derivativesUsingSequence(typename SequenceGenerator<InputsNumber>::type(),
                        sumOfDOutputs, input1, inputs...);
                }
                virtual ostream & toString(ostream & os) const { os << "cwiseProd"; return os; }

            private:
                template <int ...S>
                inline void derivativesUsingSequence(
                    Sequence<S...> seq,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input1,
                    OriginalAndDerivativeExpression<Ts> ... inputs) const {
                    // tuple of original input expressions
                    auto inputsTuple = std::make_tuple(input1.first, inputs.first ...);
                    std::tie(input1.second, inputs.second ...) =
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

                    using TraitsType = CWiseProductionTraits<typename _ForEachInputToComposeDInputI<S, InputIdx, InputExprsTupleT>::Type...>;
                    return ComposeExpression(TraitsType(), 
                        _ForEachInputToComposeDInputI<S, InputIdx, InputExprsTupleT>::expression(sumOfDOutputs, inputs)...).eval();
                }

            };





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

            template <class T, class ...Ts>
            struct ProductionResult {
                using ScalarType = scalar_type<ResultType<T>>;
                static_assert(Sequence<std::is_same<scalar_type<ResultType<Ts>>, ScalarType>::value...>::All,
                    "all scalar types of Ts MUST be the same!");
                using type =
                    decltype(ProdAll(std::declval<ResultType<T>>(), std::declval<ResultType<Ts>>() ...));
            };
            template <class T, class ...Ts>
            using ProductionResultType = typename ProductionResult<T, Ts...>::type;

                       
            // genral product traits
            template <class T, class ...Ts>
            struct ProductionTraits : public OpTraitsBase<ProductionResultType<T, Ts...>, T, Ts...> {
                inline OutputType value(ResultType<T> input1, ResultType<Ts> ... inputs) const {
                    return ProdAll(input1, inputs...);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input1,
                    OriginalAndDerivativeExpression<Ts> ... inputs) const {
                    derivativesUsingSequence(typename SequenceGenerator<InputsNumber>::type(),
                        sumOfDOutputs, input1, inputs...);
                }
                virtual ostream & toString(ostream & os) const { os << "prod"; return os; }

            private:
                template <int ...S>
                inline void derivativesUsingSequence(
                    Sequence<S...> seq,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<T> input1,
                    OriginalAndDerivativeExpression<Ts> ... inputs) const {
                    auto inputsTuple = std::make_tuple(input1.first, inputs.first ...);
                    std::tie(input1.second, inputs.second ...) = 
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

                    using TraitsType = ProductionTraits<typename 
                        _ForEachInputToComposeDInputI<S, InputIdx, InputExprsTupleT>::Type ...>;
                    auto dinput = ComposeExpression(TraitsType(),
                        _ForEachInputToComposeDInputI<S, InputIdx, InputExprsTupleT>::expression(sumOfDOutputs, inputs)...).eval();
                    return expression_assign_to<DerivativeType<InputType<InputIdx>>>(dinput);
                }
            };

        }

        template <class T, class ...Ts>
        inline Expression<CWiseProductionResultType<T, Ts...>> cwise_product(const Expression<T> & a, const Expression<Ts> & ...inputs) {
            return ComposeExpression(CWiseProductionTraits<T, Ts...>(), a, inputs...);
        }

        template <class T, class ...Ts>
        inline Expression<ProductionResultType<T, Ts...>> general_product(const Expression<T> & a, const Expression<Ts> & ...inputs) {
            return ComposeExpression(ProductionTraits<T, Ts...>(), a, inputs...);
        }

        template <class T1, class T2>
        inline Expression<ProductionResultType<T1, T2>> operator * (const Expression<T1> & a, const Expression<T2> & b) {
            return general_product(a, b);
        }




        







        // the uncyclic directional expression graph
        class ExpressionGraph {
        public:
            struct EHandleComp {
                inline bool operator()(EHandle a, EHandle b){
                    return a.id < b.id;
                }
            };

            struct IsForwardConnectionPred {
                inline IsForwardConnectionPred(const ExpressionGraph & g) :graph(g){}
                inline bool operator()(CHandle ch) const { return graph.isForwardConnection(ch); }
                const ExpressionGraph & graph;
            };

            struct IsBackwardConnectionPred {
                inline IsBackwardConnectionPred(const ExpressionGraph & g) :graph(g){}
                inline bool operator()(CHandle ch) const { return graph.isBackwardConnection(ch); }
                const ExpressionGraph & graph;
            };

        public:

            // the internal graph structure
            inline const GraphType & graph() const { return _g; }

            // input handles
            std::vector<EHandle> inputs(EHandle h) const;

            // get expression op
            inline const Op & op(EHandle h) const { return *_g.data(h).get(); }

            // get connected expression handle
            inline EHandle connected(CHandle ch) const { return _g.topo(ch).to(); }

            // check whether this is a forward connection, old -> new
            bool isForwardConnection(CHandle h) const;
            // check whether this is a backward connection, new -> old
            inline bool isBackwardConnection(CHandle h) const { return !isForwardConnection(h); }

            // get all forward connections to retrieve outputed expressions
            inline ConstConditionalContainerWrapper<HandleArray<HalfTopo>, IsForwardConnectionPred> 
                forwardConnections(EHandle h) const  {
                return MakeConditionalContainer(&(_g.topo(h).halfedges), IsForwardConnectionPred(*this));
            }
            // get all backward connections to retrieve inputed expressions
            inline ConstConditionalContainerWrapper<HandleArray<HalfTopo>, IsBackwardConnectionPred>
                backwardConnections(EHandle h) const  {
                return MakeConditionalContainer(&(_g.topo(h).halfedges), IsBackwardConnectionPred(*this));
            }

            // add new node
            EHandle addNode(std::shared_ptr<Op> op, const std::vector<EHandle>& inputs = std::vector<EHandle>());

            // add a constant value expression 
            template <class T>
            inline Expression<storage_type<T>> addConst(const T & d, const std::string & nm = "") {
                return as<storage_type<T>>(addNode(std::make_shared<OpWithConstant<storage_type<T>>>(d, nm)));
            }

            // add a reference expression
            template <class T>
            inline Expression<const T&> addRef(T & p, const std::string & nm = "") {
                return as<const T&>(addNode(std::make_shared<OpWithReference<T>>(p, nm)));
            }

            template <class T>
            inline Expression<T> as(EHandle h) {
                return Expression<T>(_g.data(h).get());
            }

            template <class T>
            inline DerivativeExpression<T> asDerivative(EHandle h) {
                return DerivativeExpression<T>(_g.data(h).get());
            }

            // evaluate the expression based on current values of vars
            void forwardPropagateExecute(EHandle result, const std::set<EHandle, EHandleComp> & vars) const;

            // create derivative graph
            // returns the derivative nodes of vars
            std::vector<EHandle> backPropagateGradient(EHandle cost, const std::vector<EHandle> & vars);

            // get expression string
            std::ostream & toString(std::ostream & os, EHandle h) const;

        private:
            GraphType _g;
        };




        // the expression wrapper
        template <class T>
        struct Expression {
        public:
            using Type = T;
            using StorageType = storage_type<T>;
            inline explicit Expression(Op * op = nullptr) : _op(op){}

            inline bool isValid() const { return _op && _op->graph && _op->self.isValid(); }

            inline EHandle handle() const { return _op ? _op->self : EHandle(); }
            inline ExpressionGraph * g() const { return _op ? _op->graph : nullptr; }

            // the current value in expression data
            inline bool hasValue() const { return _op->hasValue<T>(); }
            inline T value() const { return _op->value<T>(); }
            inline bool hasCache() const { return _op->hasCache<T>(); }
            inline const StorageType & cache() const { return _op->cache<T>(); }

            inline ResultType<T> result() const { return TraitsAboutOp<T>::Result(*_op); }

            // execute the expression graph to get value updated
            template <class ... VarTs>
            inline ResultType<T> execute(const Expression<VarTs> &... vars) const {
                assert(isValid() && "Invalid expression!");
                executeUsingSequence(typename SequenceGenerator<sizeof...(VarTs)>::type(), vars...);
                return result();
            }

            template <class ... VarTs>
            inline std::tuple<DerivativeExpression<VarTs>...> derivatives(const Expression<VarTs> &... vars){
                assert(traits<T>::is_scalar && "The cost function must has a scalar value!");
                assert(isValid() && "Invalid expression!");
                return backPropagateGradientUsingSequence(typename SequenceGenerator<sizeof...(VarTs)>::type(), vars...);
            }

            // expression cast
            template <class K>
            inline Expression<K> castTo() const { return expression_cast<K>(*this); }
            // expression eval
            inline Expression<storage_type<T>> eval() const { return expression_eval(*this); }
            // expression assign
            /*inline Expression & operator = (const Expression & e) { _op = e._op; return *this; }
            template <class K>
            inline Expression & operator = (const Expression<K> & e) { 
                return *this = ComposeExpression(ExpressionAssignTraits<T, K>(), e);
            }*/

            // sum 
            inline Expression<scalar_type<T>> sum() const { return sum_elements(*this); }

        private:
            template <class ... VarTs, int ...S>
            inline void executeUsingSequence(Sequence<S...>, const Expression<VarTs> &... vars) const {
                _op->graph->forwardPropagateExecute(_op->self, { vars.handle()... });
            }

            template <class ... VarTs, int ...S>
            inline std::tuple<DerivativeExpression<VarTs>...> backPropagateGradientUsingSequence(Sequence<S...>,
                const Expression<VarTs> &... vars){
                auto derivHandles = _op->graph->backPropagateGradient(_op->self, std::vector<EHandle>{vars.handle()...});
                return std::make_tuple(_op->graph->asDerivative<VarTs>(derivHandles[S])...);
            }

        private:
            Op * _op;
        };

        



    }
}



namespace std {

    template <typename T>
    inline ostream & operator << (ostream & s, const panoramix::core::Expression<T> & e) {
        return e.g()->toString(s, e.handle());
    }

}






 
#endif