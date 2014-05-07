#ifndef PANORAMIX_CORE_EXPRESSION_HPP
#define PANORAMIX_CORE_EXPRESSION_HPP

#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>
#include <numeric>
#include <memory>
#include <set>
#include <initializer_list>

#include "template_utilities.hpp"
#include "mesh.hpp"

#include "data_traits.hpp"


namespace panoramix {
    namespace deriv {

        using std::ostream;

        struct Op;
        using GraphType = Mesh<std::shared_ptr<Op>>;
        using EHandle = GraphType::VertHandle;
        using CHandle = GraphType::HalfHandle;

        class ExpressionGraph;

        template <class T> struct OpWithValue;
        template <class T> struct OpWithCache;
        template <class T> struct Expression;

        template <class T> using DerivativeType = DataStorageType<T>;
        template <class T> using DerivativeExpression = Expression<DataStorageType<T>>;

        // op
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

            // compute the sum of input doutputs
            // requires that all the doutputs are of the same type, 
            // and the type should be deductible during compilation
            virtual EHandle homomorphicSum(const std::vector<EHandle> & doutputs) const = 0;

            // make input derivative expressions based on inputs and output derivative expressions
            // the input derivatives should be calculated based on 
            //  1. the sum of outputDerivs 
            //  2. the inputs
            // returns
            //  1. the sum of outputDerivs
            //  2. derivative Ops corresponding to inputs
            // dinputs size MUST be same with inputs number
            virtual std::vector<EHandle> backPropagateGradient(EHandle sumOfDOutputs) const = 0;

            // make a new Op which outputs a scalar one
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
            const DataStorageType<T> & cache() const {
                return static_cast<const OpWithCache<T> *>(this)->cache;
            }

        };





        //////////////////////////////////////////////////////////////////////////////
        // DataTraits about op
        struct ResultRetrievedByValueTag {};
        struct ResultRetrievedByCacheTag {};
        template <class T>
        struct ResultTag {
            using type =
            std::conditional_t < std::is_lvalue_reference<T>::value, ResultRetrievedByValueTag,
            std::conditional_t < DataTraits<T>::shouldBeCached, ResultRetrievedByCacheTag,
            ResultRetrievedByValueTag >> ;
        };

        template <class T, class ResultTag = typename ResultTag<T>::type>
        struct TraitsAboutOp {};

        template <class T> // by value
        struct TraitsAboutOp<T, ResultRetrievedByValueTag> {
            using OpType = OpWithValue<T>;
            using ResultType = T;
            static T Result(Op const & op) { // use value<T>() directly to fetch result
                if (!op.hasValue<T>()){
                    throw std::runtime_error(std::string("it is determined that ") + typeid(T).name() +
                        " is retrieved by value(), but current op doesn't has the correct value!");
                }
                return op.value<T>();
            }
        };

        template <class T>
        struct TraitsAboutOp<T, ResultRetrievedByCacheTag> { // by cache
            using OpType = OpWithCache<T>;
            using ResultType = const DataStorageType<T> &;
            static ResultType Result(Op const & op) {  // use reference to DataStorageType to fetch result
                if (!op.hasCache<T>()){
                    throw std::runtime_error(std::string("it is determined that ") + typeid(T).name() +
                        " is retrieved by cache(), but current op doesn't has the correct cache!");
                }
                return op.cache<T>();
            }
        };


        // smartly determines whether cache should be used for this type
        template <class T>
        using OpBaseType = typename TraitsAboutOp<T>::OpType;

        // smartly determines the result return type
        template <class T>
        using ResultType = typename TraitsAboutOp<T>::ResultType;






        /////////////////////////////////////////////////////////////////////////////////////
        // op with typed value
        template <class T>
        struct OpWithValue : public Op {

            // compute value
            virtual T value() const = 0;

            // do nothing
            virtual void forwardPropagateExecute() {}

            // make a scalar one to initialize automatic derivation
            virtual std::shared_ptr<Op> one() const override {
                return std::make_shared<OpWithConstant<DerivativeType<T>>>(DataTraits<DerivativeType<T>>::one());
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
            using ST = DataStorageType<T>;
            static_assert(std::is_assignable<ST&, T>::value, "T cannot be assigned to ST&!");
            inline OpWithCache(const ST & c = ST(), const std::string nm = "")
                : cache(c), givenName(nm)
            {}

            virtual std::ostream & toString(std::ostream & os) const override { 
                os << givenName; 
                return os; 
            }

            // compute value and store in cache
            virtual void forwardPropagateExecute() { cache = DataTraits<T>::eval(value()); }

            ST cache;
            std::string givenName;
        };

        // constant
        template <class T>
        struct OpWithConstant : public OpWithCache<T> {
            static_assert(IsStorageType<T>::value, "T must be a storage type!");
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
        inline OriginalAndDerivativeExpression<T>
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
                        graph->asDerived<OutputType>(sumOfDOutputs), // DOutput
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

        // expression cast
        namespace  {
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
            struct ExpressionCastTraits : public OpTraitsBase<DataStorageType<To>, From> {
                inline DataStorageType<To> value(ResultType<From> from) const {
                    return DataTraits<To>::castFromWithScalarConversion<ResultType<From>>(from);
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

            // eval
            template <class T>
            struct ExpressionEvalTraits : public OpTraitsBase<DataStorageType<ResultType<T>>, T> {
                inline OutputType value(ResultType<T> from) const {
                    return DataTraits<ResultType<T>>::eval(std::move(from));
                }
                static_assert(std::is_same<DerivativeExpression<OutputType>,
                    DerivativeExpression<T >> ::value, "should be the same!");
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
        std::enable_if_t<std::is_same<To, From>::value, Expression<To>> 
            expressionAssign(const Expression<From> & from) {
            return from;
        }

        template <class To, class From>
        std::enable_if_t<!std::is_same<To, From>::value, Expression<To>> 
            expressionAssign(const Expression<From> & from) {
            return ComposeExpression(ExpressionAssignTraits<To, From>(), from);
        }

        template <class To, class From>
        std::enable_if_t<std::is_same<DataStorageType<To>, From>::value, 
            Expression<DataStorageType<To>>> 
            expressionCast(const Expression<From> & from) {
            return from;
        }

        template <class To, class From>
        std::enable_if_t<!std::is_same<DataStorageType<To>, From>::value, Expression<DataStorageType<To>>>
            expressionCast(const Expression<From> & from) {
            return ComposeExpression(ExpressionCastTraits<To, From>(), from);
        }

        template <class T>
        std::enable_if_t<std::is_same<DataStorageType<T>, T>::value, 
            Expression<DataStorageType<T>>> 
            expressionEval(const Expression<T> & from) {
            return from;
        }

        template <class T>
        std::enable_if_t<!(std::is_same<DataStorageType<T>, T>::value), 
            Expression<DataStorageType<T>>>
            expressionEval(const Expression<T> & from) {
            return ComposeExpression(ExpressionEvalTraits<T>(), from);
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

            // the  graph structure
            inline const GraphType & graph() const { return _g; }

            // reserve expressions number
            void reserve(size_t sz);

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
            inline Expression<DataStorageType<T>> addConst(const T & d, const std::string & nm = "") {
                return as<DataStorageType<T>>(addNode(std::make_shared<OpWithConstant<DataStorageType<T>>>(d, nm)));
            }

            // add a reference expression
            template <class T>
            inline Expression<const T&> addRef(T & p, const std::string & nm = "") {
                return as<const T&>(addNode(std::make_shared<OpWithReference<T>>(p, nm)));
            }

            // make an expression from a handle
            template <class T>
            inline Expression<T> as(EHandle h) {
                return Expression<T>(_g.data(h).get());
            }

            // make a derivative expression from a handle
            template <class T>
            inline DerivativeExpression<T> asDerived(EHandle h) {
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





        // the expression class
        // a wrapper of an OpWithValue
        template <class T>
        struct Expression {
        public:
            using Type = T;
            inline explicit Expression(Op * op = nullptr) : _op(op){}

            inline bool isValid() const { return _op && _op->graph && _op->self.isValid(); }

            inline EHandle handle() const { return _op ? _op->self : EHandle(); }
            inline ExpressionGraph * g() const { return _op ? _op->graph : nullptr; }

            // the current value in expression data
            inline ResultType<T> result() const { return TraitsAboutOp<T>::Result(*_op); }

            // execute the expression graph to get value updated
            template <class ... VarTs>
            inline ResultType<T> execute(const Expression<VarTs> &... vars) const {
                assert(isValid() && "Invalid expression!");
                executeUsingSequence(typename SequenceGenerator<sizeof...(VarTs)>::type(), vars...);
                return result();
            }

            // execute the expression graph to get value updated
            template <class ExpressionIteratorT>
            inline ResultType<T> execute(ExpressionIteratorT begin, ExpressionIteratorT end) const {
                assert(isValid() && "Invalid expression!");
                std::set<EHandle, ExpressionGraph::EHandleComp> dependencies;
                for (auto i = begin; i != end; ++i)
                    dependencies.insert(i->handle());
                _op->graph->forwardPropagateExecute(_op->self, dependencies);
                return result();
            }

            // compute the derivative of current expression
            template <class ... VarTs>
            inline std::tuple<DerivativeExpression<VarTs>...> derivatives(const Expression<VarTs> &... vars) const {
                assert(IsScalarType<T>::value && "The cost function must has a scalar value!");
                assert(isValid() && "Invalid expression!");
                return backPropagateGradientUsingSequence(typename SequenceGenerator<sizeof...(VarTs)>::type(), vars...);
            }

            // compute the derivative of current expression
            template <class ExpressionIteratorT, class ExpressionOutIteratorT>
            inline void derivatives(ExpressionIteratorT varsBegin, ExpressionIteratorT varsEnd, ExpressionOutIteratorT derivsBegin) const {
                assert(IsScalarType<T>::value && "The cost function must has a scalar value!");
                assert(isValid() && "Invalid expression!");

                using InputsType = typename std::iterator_traits<ExpressionIteratorT>::value_type::Type;
                using OutputsType = DerivativeType<InputsType>;

                std::vector<EHandle> varHandles;
                for (auto i = varsBegin; i != varsEnd; ++i)
                    varHandles.push_back(i->handle());
                auto derivHandles = _op->graph->backPropagateGradient(_op->self, varHandles);

                std::transform(derivHandles.begin(), derivHandles.end(), derivsBegin, 
                    [this](EHandle h){
                    return _op->graph->as<OutputsType>(h); 
                });
            }

            // expression assign
            template <class K>
            inline Expression<K> assign() const { return expressionAssign<K>(*this); }

            // expression cast
            template <class K>
            inline Expression<DataStorageType<K>> cast() const { return expressionCast<K>(*this); }

            // expression eval
            inline Expression<DataStorageType<T>> eval() const {
                return expressionEval<T>(*this); 
            }


            void printType() const { std::cout << typeid(T).name() << std::endl; }


            inline Expression<DataScalarType<T>> sum() const { return sumElements(*this); }

        private:
            template <class ... VarTs, int ...S>
            inline void executeUsingSequence(Sequence<S...>, const Expression<VarTs> &... vars) const {
                _op->graph->forwardPropagateExecute(_op->self, { vars.handle()... });
            }

            template <class ... VarTs, int ...S>
            inline std::tuple<DerivativeExpression<VarTs>...> backPropagateGradientUsingSequence(Sequence<S...>,
                const Expression<VarTs> &... vars) const {
                auto derivHandles = _op->graph->backPropagateGradient(_op->self, std::vector<EHandle>{vars.handle()...});
                return std::make_tuple(_op->graph->asDerived<VarTs>(derivHandles[S])...);
            }

        private:
            Op * _op;
        };

        



    }
}



namespace std {

    template <typename T>
    inline ostream & operator << (ostream & s, const panoramix::deriv::Expression<T> & e) {
        return e.g()->toString(s, e.handle());
    }

}






 
#endif