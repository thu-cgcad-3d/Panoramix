#ifndef PANORAMIX_CORE_EXPRESSION_HPP
#define PANORAMIX_CORE_EXPRESSION_HPP

#include <iostream>
#include <utility>
#include <algorithm>
#include <numeric>
#include <initializer_list>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "basic_types.hpp"
#include "template_utilities.hpp"
#include "mesh.hpp"
#include "data_traits.hpp"
#include "simple_data_traits.hpp"

namespace panoramix {
    namespace core {


        class ExpressionGraph;
        struct Op;

        using GraphType = Mesh<std::shared_ptr<Op>>;
        using EHandle = GraphType::VertHandle;
        using CHandle = GraphType::HalfHandle;

        // the expression wrapper
        template <class DataT>
        struct Expression {
            static_assert(IsActive<DataT>::value, "type is not currently supported!");

        public:
            using Type = DataT;
            inline Expression(EHandle h = EHandle(), ExpressionGraph * g = nullptr) : _h(h), _g(g){}

            inline EHandle handle() const { return _h; }
            inline ExpressionGraph * g() const { return _g; }

            // the current value in expression data
            inline const DataT & value() const { return _g->op(_h).as<DataT>().value(); }

            // evaluate the expression and return the newest value
            template <class ... VarTs>
            inline const DataT & eval(const Expression<VarTs> &... vars) const { 
                assert(_h.isValid() && "Invalid expression handle!");
                return evalUsingSequence(typename SequenceGenerator<sizeof...(VarTs)>::type(), vars...);
            }

            template <class ... VarTs>
            inline std::tuple<Expression<VarTs>...> derivatives(const Expression<VarTs> &... vars){
                assert(GetSize(value()) == 1 && "The cost function must has a scalar value!");
                assert(_h.isValid() && "Invalid expression handle!");
                return derivativesUsingSequence(typename SequenceGenerator<sizeof...(VarTs)>::type(), vars...);
            }

        private:
            template <class ... VarTs, int ...S>
            inline const DataT & evalUsingSequence(Sequence<S...>, const Expression<VarTs> &... vars) const {
                _g->evaluate(_h, {vars.handle()...});
                return value();
            }

            template <class ... VarTs, int ...S>
            inline std::tuple<Expression<VarTs>...> derivativesUsingSequence(Sequence<S...>, 
                const Expression<VarTs> &... vars){
                auto derivHandles = _g->createDerivatives(_h, std::vector<EHandle>{vars.handle()...});
                return std::make_tuple(Expression<VarTs>(derivHandles[S], _g)...);
            }            

        private:
            EHandle _h;
            ExpressionGraph * _g;
        };


        template <class DataT>
        struct OpWithAValue;
        
        // operator base
        struct Op {

            // the graph
            ExpressionGraph * graph;

            // the handle of self node
            EHandle self;

            // stream to string
            virtual std::ostream & toString(std::ostream & os) const = 0;

            // evaluate and get result
            virtual void eval(ExpressionGraph const * const g, std::vector<EHandle> && inputs) {}

            // make input derivative expressions based on inputs and output derivative expressions
            // returns derivative Ops corresponding to inputs
            virtual std::vector<EHandle> makeInputDerivatives(ExpressionGraph * const g,
                EHandle self,
                std::vector<EHandle> && inputs,
                EHandle outputDerivsSum) const {
                return std::vector<EHandle>();
            }

            // make a new Op which accepts the current Op as input, 
            // reserve the input Op's shape, and its function is to fill the shape with ones
            virtual std::shared_ptr<Op> makeOnes() const = 0;

            // make a new Op which accepts multiple inputs all with the same type the current Op, 
            // its function is to sum all the inputs
            virtual std::shared_ptr<Op> makePlus() const = 0; // support > 2 args, all inputs must share the same type

            // convert to op with a value to retrieve a value
            template <class DataT>
            const OpWithAValue<DataT> & as() const { 
                return *static_cast<const OpWithAValue<DataT> *>(this); 
            }

            // check whether this op has the value type
            template <class DataT>
            bool has() const { 
                return dynamic_cast<const OpWithAValue<DataT> *>(this) != nullptr; 
            }

        };

        // base op for expressions
        template <class DataT>
        struct OpWithAValue : public Op {
            virtual const DataT & value() const = 0;
            virtual std::shared_ptr<Op> makeOnes() const override { return std::make_shared<SetConstantOp<DataT>>(value(), ElementType<DataT>(1)); }
            virtual std::shared_ptr<Op> makePlus() const override { return std::make_shared<ElementWiseSumOp<DataT>>(value()); }
        };

        // using a value cache to store data
        template <class DataT>
        struct OpWithACache : public OpWithAValue<DataT> {
            inline OpWithACache(const DataT & d, const std::string nm = "") 
                : cache(d), givenName(nm)
            {}
            virtual std::ostream & toString(std::ostream & os) const { os << "{" << givenName << "}"; return os; }
            virtual const DataT & value() const override { return cache; }            
            DataT cache;
            std::string givenName;
        };

        // using a reference type to reference data
        template <class DataT>
        struct OpWithAReference : public OpWithAValue<DataT> {
            inline OpWithAReference(DataT & m, const std::string & nm = "") 
                : ref(m), givenName(nm)
            {}
            virtual std::ostream & toString(std::ostream & os) const { os << "{" << givenName << "&}"; return os; }
            virtual const DataT & value() const override { return ref; }
            DataT & ref;
            std::string givenName;
        };

        // using delayed computation
        template <class DataT, class DirectOutputT>
        struct OpDelayed : public OpWithACache<DataT> {
            //virtual DirectOutputT directValue(const std::vector<EHandle>)
        };


        // basic utility Ops

        // use the shape of the input value and fill it with a scalar constant
        template <class DataT>
        struct SetConstantOp : public OpWithACache<DataT> {
            using ScalarType = ElementType<DataT>;

            inline SetConstantOp(const DataT & init, const ScalarType & ss) : OpWithACache<DataT>(init), s(ss){}

            virtual std::ostream & toString(std::ostream & os) const { os << "SetConstant[" << s << "]"; return os; }

            virtual void eval(ExpressionGraph const * const g, std::vector<EHandle> && inputs) {
                assert(inputs.size() == 1);
                assert(g->op(inputs.front()).has<DataT>());
                cache = g->op(inputs.front()).as<DataT>().value(); // get shape information
                Fill(cache, s); // Fill with scalar
            }

            virtual std::vector<EHandle> makeInputDerivatives(ExpressionGraph * const g,
                EHandle self,
                std::vector<EHandle> && inputs,
                EHandle outputDerivsSum) const {
                assert(inputs.size() == 1);
                assert(g->op(outputDerivsSum).has<DataT>());
                // create a SetConstant(~, 0) for input since the values in input do not affect the output at all
                // since SetConstant is used as the source for derivative graph construction,
                // this indicates that the constructed derivative graph is connected with the original expression graph
                return std::vector<EHandle>{g->addNode(std::make_shared<SetConstantOp<DataT>>(value(), ScalarType(0)), { outputDerivsSum })};
            }

            ScalarType s;
        };

        template <class DataT>
        Expression<DataT> SetConstant(const Expression<DataT> & e, const ElementType<DataT> & s) {
            return Expression<DataT>(e.g()->addNode(std::make_shared<SetConstantOp<DataT>>(e.value(), s), { e.handle() }), e.g());
        }


        // elementwise plus
        template <class DataT>
        struct ElementWiseSumOp : public OpWithACache<DataT> {
            inline ElementWiseSumOp(const DataT & init) : OpWithACache<DataT>(init){}

            virtual std::ostream & toString(std::ostream & os) const { os << "ElementWiseSum"; return os; }

            // evaluate and get result
            virtual void eval(ExpressionGraph const * const g, std::vector<EHandle> && inputs) {
                Fill<DataT>(cache, 0);
                for (auto h : inputs){
                    cache += g->op(h).as<DataT>().value();
                }
            }

            // make input derivative expressions based on inputs and output derivative expressions
            // returns derivative Ops corresponding to inputs
            virtual std::vector<EHandle> makeInputDerivatives(ExpressionGraph * const g,
                EHandle self,
                std::vector<EHandle> && inputs,
                EHandle outputDerivsSum) const {
                std::vector<EHandle> inputDerivs(inputs.size(), outputDerivsSum);
                return inputDerivs;
            }
        };

        template <class DataT1, class ... DataTs>
        Expression<DataT1> ElementWiseSum(const Expression<DataT1> & a, const Expression<DataTs> & ... others) {
            static_assert(Sequence<std::is_same<DataTs, DataT1>::value...>::All,
                "all input expressions must share the same type!");
            ExpressionGraph * g = a.g();
            // compute current sum value
            DataT1 vals[] = { a.value(), others.value()... };
            DataT1 sumVal;
            Fill<DataT1>(sumVal, 0);
            for (DataT1 & v : vals){
                sumVal += v;
            }
            
            EHandle h = g->addNode(std::make_shared<ElementWiseSumOp<DataT1>>(sumVal),
                {a.handle(), others.handle()...});
            return Expression<DataT1>(h, g);
        }

        
        // elementwise multiplication of multiple inputs
        template <class DataT>
        struct ElementWiseProductOp : public OpWithACache<DataT> {
            inline ElementWiseProductOp(const DataT & init) : OpWithACache<DataT>(init){}

            virtual std::ostream & toString(std::ostream & os) const { os << "ElementWiseProduct"; return os; }

            // evaluate and get result
            virtual void eval(ExpressionGraph const * const g, std::vector<EHandle> && inputs) {
                Fill<DataT>(cache, 1);
                for (auto h : inputs){
                    assert(g->op(h).has<DataT>());
                    ElementWiseMult<DataT>(cache, g->op(h).as<DataT>().value(), cache);
                }
            }

            // make input derivative expressions based on inputs and output derivative expressions
            // returns derivative Ops corresponding to inputs
            virtual std::vector<EHandle> makeInputDerivatives(ExpressionGraph * const g,
                EHandle self,
                std::vector<EHandle> && inputs,
                EHandle outputDerivsSum) const {

                std::vector<EHandle> inputDerivs(inputs.size());
                for (int i = 0; i < inputs.size(); i++){
                    std::vector<EHandle> toMults(1, outputDerivsSum);
                    toMults.reserve(inputs.size());
                    for (int j = 0; j < inputs.size(); j++){
                        if (j == i)
                            continue;
                        toMults.push_back(inputs[j]);
                    }
                    inputDerivs[i] = toMults.size() == 1 ? toMults.front() :
                        g->addNode(std::make_shared<ElementWiseProductOp<DataT>>(value()), toMults);
                }

                return inputDerivs;
            }
        };

        template <class DataT1, class ... DataTs>
        Expression<DataT1> ElementWiseProduct(const Expression<DataT1> & a, const Expression<DataTs> & ... others) {
            static_assert(Sequence<std::is_same<DataTs, DataT1>::value...>::All,
                "all input expressions must share the same type!");
            ExpressionGraph * g = a.g();
            // compute current product value
            DataT1 vals[] = { a.value(), others.value()... };
            DataT1 prodVal;
            Fill<DataT1>(prodVal, 1);
            for (DataT1 & v : vals){
                ElementWiseMult<DataT1>(prodVal, v, prodVal);
            }

            EHandle h = g->addNode(std::make_shared<ElementWiseProductOp<DataT1>>(prodVal), 
                { a.handle(), others.handle()... });
            return Expression<DataT1>(h, g);
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
            template <class DataT>
            inline Expression<EvaluatedType<DataT>> addConst(const DataT & d, const std::string & nm = "") {
                return Expression<EvaluatedType<DataT>>(
                    addNode(std::make_shared<OpWithACache<EvaluatedType<DataT>>>(d, nm)), 
                    this);
            }

            // add a reference expression
            template <class DataT>
            inline Expression<EvaluatedType<DataT>> addRef(DataT & p, const std::string & nm = "") {
                return Expression<EvaluatedType<DataT>>(
                    addNode(std::make_shared<OpWithAReference<EvaluatedType<DataT>>>(p, nm)),
                    this);
            }

            // evaluate the expression based on current values of vars
            void evaluate(EHandle result, const std::set<EHandle, EHandleComp> & vars) const;

            // create derivative graph
            // returns the derivative nodes of vars
            std::vector<EHandle> createDerivatives(EHandle cost, const std::vector<EHandle> & vars); 

            // get expression string
            std::ostream & toString(std::ostream & os, EHandle h) const;

        private:
            GraphType _g;
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