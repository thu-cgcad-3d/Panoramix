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

#include "../core/template_utilities.hpp"
#include "../core/misc.hpp"
#include "../core/mesh.hpp"

#include "data_traits_definitions.hpp"

namespace panoramix {
    namespace deriv {

        using std::ostream;
        using core::Sequence;
        using core::SequenceGenerator;

        struct Op;
        using GraphType = core::Mesh<std::shared_ptr<Op>>;
        using EHandle = GraphType::VertHandle;
        using CHandle = GraphType::HalfHandle;

        class ExpressionGraph;

        template <class T> struct OpWithValue;
        template <class T> struct OpWithCache;
        template <class T> struct Expression;

        template <class T> struct IsExpression : public std::false_type {};
        template <class T> struct IsExpression<Expression<T>> : public std::true_type {};

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
            virtual void forwardPropagateExecute() { cache = common::Eval(value()); }

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
        

        




        // expression assign
        template <class To, class From>
        inline std::enable_if_t<std::is_same<To, From>::value, Expression<To>> 
            expressionAssign(const Expression<From> & from) {
            return from;
        }

        template <class To, class From>
        inline std::enable_if_t<!std::is_same<To, From>::value, Expression<To>>
            expressionAssign(const Expression<From> & from) {
            struct AssignerOp : public OpBaseType<To> {
                virtual std::ostream & toString(std::ostream & os) const {
                    os << "assign";
                    return os;
                }
                virtual To value() const override {
                    std::vector<EHandle> inp = graph->inputs(self);
                    assert(inp.size() == 1 && "wrong inputs number");
                    return TraitsAboutOp<From>::Result(graph->op(inp.front()));
                }
                virtual std::vector<EHandle> backPropagateGradient(EHandle sumOfDOutputs) const {
                    return std::vector<EHandle>(1, 
                        expressionAssign<DerivativeType<From>, DerivativeType<To>>(
                        graph->as<DerivativeType<To>>(sumOfDOutputs)).handle());
                }
            };
            return from.g()->as<To>(from.g()->addNode(std::make_shared<AssignerOp>(), {from.handle()}));
        }

        // expression cast
        template <class To, class From>
        inline std::enable_if_t<std::is_same<To, From>::value,
            Expression<To>> 
            expressionCast(const Expression<From> & from) {
            return from;
        }

        template <class To, class From>
        inline std::enable_if_t<!std::is_same<To, From>::value, Expression<To>> 
            expressionCast(const Expression<From> & from) {
            struct CasterOp : public OpBaseType<To> {
                virtual std::ostream & toString(std::ostream & os) const {
                    os << "cast";
                    return os;
                }
                virtual To value() const override {
                    std::vector<EHandle> inp = graph->inputs(self);
                    assert(inp.size() == 1 && "wrong inputs number");
                    To to;
                    common::Cast(TraitsAboutOp<From>::Result(graph->op(inp.front())), to);
                    return to;
                    //return DataTraits<To>::castFromWithScalarConversion(TraitsAboutOp<From>::Result(graph->op(inp.front())));
                }
                virtual std::vector<EHandle> backPropagateGradient(EHandle sumOfDOutputs) const {
                    return std::vector<EHandle>(1,
                        expressionCast<DerivativeType<From>, DerivativeType<To>>(
                        graph->as<DerivativeType<To>>(sumOfDOutputs)).handle());
                }
            };
            return from.g()->as<To>(from.g()->addNode(std::make_shared<CasterOp>(), { from.handle() }));
        }

        // expression eval
        template <class T>
        inline std::enable_if_t<std::is_same<DataStorageType<T>, T>::value,
            Expression<DataStorageType<T>>> 
            expressionEval(const Expression<T> & from) {
            return from;
        }

        template <class T>
        inline std::enable_if_t<!(std::is_same<DataStorageType<T>, T>::value), 
            Expression<DataStorageType<T>>>
            expressionEval(const Expression<T> & from) {
            struct EvalOp : public OpBaseType<DataStorageType<T>> {
                virtual std::ostream & toString(std::ostream & os) const {
                    os << "eval";
                    return os;
                }
                virtual DataStorageType<T> value() const override {
                    std::vector<EHandle> inp = graph->inputs(self);
                    assert(inp.size() == 1 && "wrong inputs number");
                    //return DataTraits<T>::eval(TraitsAboutOp<T>::Result(graph->op(inp.front())));
                    return common::Eval(TraitsAboutOp<T>::Result(graph->op(inp.front())));
                }
                virtual std::vector<EHandle> backPropagateGradient(EHandle sumOfDOutputs) const {
                    return std::vector<EHandle>(1, sumOfDOutputs);
                }
            };
            return from.g()->as<DataStorageType<T>>(from.g()->addNode(std::make_shared<EvalOp>(), { from.handle() }));
        }

        //// expression conditional branch
        //template <class T, class SelectorT, class ConditionT, class ExpressionTIteratorT>
        //inline Expression<T> expressionOnCondition(SelectorT && cc,
        //    const Expression<ConditionT> & cond,
        //    ExpressionTIteratorT candidatesBegin, ExpressionTIteratorT candidatesEnd) {
        //    typename InputType = typename std::iterator_traits<ExpressionTIteratorT>::value_type::Type;
        //    static_assert(std::is_same<T, InputType>::value, "invalid inputs type!");
        //    struct OnConditionOp : public OpBaseType<T> {
        //        inline explicit OnConditionOp(SelectorT && cc) : selector(std::forward<SelectorT>(cc)) {}
        //        virtual std::ostream & toString(std::ostream & os) const {
        //            os << "onCondition";
        //            return os;
        //        }
        //        virtual T value() const override {
        //            std::vector<EHandle> inp = graph->inputs(self);
        //            assert(inp.size() >= 1 && "wrong inputs number");
        //            int choice = selector(TraitsAboutOp<T>::Result(graph->op(inp.front())));
        //            assert(choice >= 0 && choice < inp.size() - 1 && "invalid choice!");
        //            return TraitsAboutOp<T>::Result(graph->op(inp[choice + 1]));
        //        }
        //        virtual std::vector<EHandle> backPropagateGradient(EHandle sumOfDOutputs) const {
        //            std::vector<EHandle> inp = graph->inputs(self);
        //            assert(inp.size() >= 1 && "wrong inputs number");
        //            int choice = selector(TraitsAboutOp<T>::Result(graph->op(inp.front())));
        //            assert(choice >= 0 && choice < inp.size() && "invalid choice!");
        //            std::vector<EHandle> dinputs(inp.size());
        //            dinputs[choice + 1] = sumOfDOutputs;
        //            return dinputs;
        //        }
        //        SelectorT selector;
        //    };
        //    std::vector<EHandle> inputs(1, cond.handle());
        //    std::transform(candidatesBegin, candidatesEnd, std::back_inserter(inputs), 
        //        [](const Expression<T> & branch){
        //        return branch.handle(); 
        //    });
        //    return cond.g()->as<T>(cond.g()
        //        ->addNode(std::make_shared<OnConditionOp>(std::forward<SelectorT>(cc)), inputs));
        //}     


       




        







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

            // invalidate all expressions
            void invalidateAll();

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
            inline core::ConstConditionalContainerWrapper<core::HandleArray<core::HalfTopo>, IsForwardConnectionPred> 
                forwardConnections(EHandle h) const  {
                return core::MakeConditionalContainer(&(_g.topo(h).halfedges), IsForwardConnectionPred(*this));
            }
            // get all backward connections to retrieve inputed expressions
            inline core::ConstConditionalContainerWrapper<core::HandleArray<core::HalfTopo>, IsBackwardConnectionPred>
                backwardConnections(EHandle h) const  {
                return core::MakeConditionalContainer(&(_g.topo(h).halfedges), IsBackwardConnectionPred(*this));
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