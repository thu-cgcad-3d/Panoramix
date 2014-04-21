#ifndef PANORAMIX_CORE_EXPRESSION_HPP
#define PANORAMIX_CORE_EXPRESSION_HPP

#include <iostream>
#include <utility>
#include <initializer_list>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "basic_types.hpp"
#include "template_utilities.hpp"
#include "mesh.hpp"
#include "data_traits.hpp"

namespace panoramix {
    namespace core {

        /*
        enum DimensionType : int {
            Dyn = 0
        };


        // shape
        template <int ... Ds>
        class Shape;

        // empty shape
        template <>
        class Shape<> {
        public:
            static const int StaticVolume = 1;
            static const int Rank = 0;
            static const int DynamicNum = 0;
            inline Shape(){}
            inline int volume() const { return 1; }
        };

        // shape with fixed dimension
        template <int D, int ... Ds>
        class Shape<D, Ds...> : private Shape<Ds...>{
        public:
            using BaseType = Shape<Ds...>;
            static const bool LengthIsDynamic = false;
            static const int StaticLength = D;
            static const int StaticVolume = BaseType::StaticVolume * D;
            static const int Rank = sizeof...(Ds)+1;
            static const int DynamicNum = BaseType::DynamicNum;

            // ignore fixed dimension in constructor
            template <class... Ints>
            inline explicit Shape(Ints ... ints)
                : BaseType(ints...){}

            // promote
            inline explicit Shape(const BaseType & s)
                : BaseType(s){}

            // base
            inline const BaseType& base() const {
                return (const BaseType&)(*this);
            }

            // first order size
            inline int length() const { return D; }

            // volume
            inline int volume() const { return ((const BaseType&)(*this)).volume() * length(); }

            // size
            template <int Ord> inline int size() const {
                return ((const BaseType&)(*this)).size<Ord - 1>();
            }
            template <> inline int size<0>() const { return length(); }

            // resize
            template <int Ord> inline void resize(int s) {
                return ((BaseType&)(*this)).resize<Ord - 1>(s);
            }
            template <> inline void resize<0>(int s) {
                static_assert(false, "unresizable on fixed dimension!");
            }

            // dynamic size
            template <int DOrd> inline int dsize() const {
                return ((const BaseType&)(*this)).dsize<DOrd>();
            }

            // dynamic resize
            template <int DOrd> inline void dresize(int ds) {
                ((BaseType&)(*this)).dresize<DOrd>(ds);
            }


            // static size
            template <int SOrd> inline int ssize() const {
                return ((const BaseType&)(*this)).ssize<SOrd - 1>();
            }
            template <> inline int ssize<0>() const {
                return D;
            }

        };

        // shape with dynamic dimension
        template <int ... Ds>
        class Shape<Dyn, Ds...> : private Shape<Ds...>{
        public:
            using BaseType = Shape<Ds...>;
            static const bool LengthIsDynamic = true;
            static const int StaticLength = Dyn;
            static const int StaticVolume = Dyn;
            static const int Rank = sizeof...(Ds)+1;
            static const int DynamicNum = BaseType::DynamicNum + 1;

            // set size for Dyn dimension
            template <class Int, class... Ints>
            inline explicit Shape(Int i, Ints ... ints)
                : BaseType(ints...), _length(i){}

            // promote
            inline Shape(int len, const Shape<Ds...> & s)
                : BaseType(s), _length(len){}



            // base
            inline const BaseType& base() const {
                return (const BaseType&)(*this);
            }

            // first order size
            inline int length() const { return _length; }

            // volume
            inline int volume() const { return ((const BaseType&)(*this)).volume() * length(); }

            // size
            template <int Ord> inline int size() const {
                return ((const BaseType&)(*this)).size<Ord - 1>();
            }
            template <> inline int size<0>() const { return _length; }

            // resize
            template <int Ord> inline void resize(int s) {
                return ((BaseType&)(*this)).resize<Ord - 1>(s);
            }
            template <> inline void resize<0>(int s) {
                _length = s;
            }

            // dynamic size
            template <int DOrd> inline int dsize() const {
                return ((const BaseType&)(*this)).dsize<DOrd - 1>();
            }
            template <> inline int dsize<0>() const {
                return _length;
            }

            // dynamic resize
            template <int DOrd> inline void dresize(int ds) {
                ((BaseType&)(*this)).dresize<DOrd - 1>(ds);
            }
            template <> inline void dresize<0>(int ds) {
                _length = ds;
            }

            // static size
            template <int SOrd> inline int ssize() const {
                return ((const BaseType&)(*this)).ssize<SOrd>();
            }


        private:
            int _length;
        };


        namespace {

            /////////////////////////
            template <int ...Ds, int ...S>
            inline std::array<int, sizeof...(Ds)> MakeArrayFromShapeUsingSequence(const Shape<Ds...> & s, 
                Sequence<S...>){
                return std::array<int, sizeof...(Ds)>{ { s.size<S>()... } };
            }

            template <int ...Ds, int ...S>
            inline std::vector<int> MakeVectorFromShapeUsingSequence(const Shape<Ds...> & s,
                Sequence<S...>){
                return std::vector<int>{ s.size<S>()... };
            }

            template <int ...Ds, int ...S>
            inline auto MakeTupleFromShapeUsingSequence(const Shape<Ds...> & s,
                Sequence<S...>) -> decltype(std::make_tuple((int)s.size<S>()...)){
                return std::make_tuple((int)s.size<S>()... );
            }

        }

        // array from shape
        template <int ...Ds>
        inline std::array<int, sizeof...(Ds)> MakeArrayFromShape(const Shape<Ds...> & s){
            return MakeArrayFromShapeUsingSequence(s, 
                typename SequenceGenerator<sizeof...(Ds)>::type());
        }

        // vector from shape
        template <int ...Ds>
        inline std::vector<int> MakeVectorFromShape(const Shape<Ds...> & s){
            return MakeVectorFromShapeUsingSequence(s,
                typename SequenceGenerator<sizeof...(Ds)>::type());
        }

        // tuple from shape
        template <int ...Ds>
        inline auto MakeTupleFromShape(const Shape<Ds...> & s) 
            -> decltype(MakeTupleFromShapeUsingSequence(s, 
            typename SequenceGenerator<sizeof...(Ds)>::type())) {
            return MakeTupleFromShapeUsingSequence(s, 
                typename SequenceGenerator<sizeof...(Ds)>::type());
        }
            
        namespace {

            ////////////////////////
            template <class SizeAndSubPair>
            inline int ComputeIndexFromSizeAndSubs(SizeAndSubPair szsb){
                return szsb.second;
            }

            template <class SizeAndSubPair, class ...SizeAndSubPairs>
            inline int ComputeIndexFromSizeAndSubs(SizeAndSubPair szsb, SizeAndSubPairs... sss){
                auto i = ComputeIndexFromSizeAndSubs(sss...) * szsb.first + szsb.second;
                return i;
            }

            template <int ...Ds, class SubsArray, int ...S>
            inline int ComputeIndexFromShapeAndSubsUsingSequence(const Shape<Ds...> & shape,
                Sequence<S...>, SubsArray subs){
                return ComputeIndexFromSizeAndSubs(std::make_pair(shape.size<S>(), std::get<S>(subs))...);
            }

            ////////////////////////
            template <class Sub>
            inline std::pair<int, Sub&> MakePairToReturn(int size, Sub& sub){
                return std::pair<int, Sub&>{size, sub};
            }

            template <class SizeAndSubPair>
            inline void ComputeSubsFromSizeAndIndex(int index, SizeAndSubPair szsb){
                szsb.second = index % szsb.first;
            }

            template <class SizeAndSubPair, class ...SizeAndSubPairs>
            inline void ComputeSubsFromSizeAndIndex(int index, SizeAndSubPair szsb, SizeAndSubPairs... sss){
                szsb.second = index % szsb.first;
                ComputeSubsFromSizeAndIndex(index / std::get<0>(szsb), sss...);
            }

            template <int ...Ds, class SubsArray, int ...S>
            inline void ComputeSubsFromShapeAndIndexUsingSequence(const Shape<Ds...> & shape,
                int index, Sequence<S...>, SubsArray & subs){
                return ComputeSubsFromSizeAndIndex(index, MakePairToReturn(shape.size<S>(), std::get<S>(subs))...);
            }

        }

        // sub2ind
        template <int ...Ds, class ...Subs>
        inline int IndexFromSubs(const Shape<Ds...> & shape, Subs... subs){
            static_assert(sizeof...(Ds) == sizeof...(Subs), "subcript length not compatible with shape rank!");
            return ComputeIndexFromShapeAndSubsUsingSequence(shape,
                typename SequenceGenerator<sizeof...(Ds)>::type(), std::tie(subs...));
        }

        template <int ...Ds>
        inline int IndexFromSubs(const Shape<Ds...> & shape, const std::array<int, sizeof...(Ds)> & subs){
            return ComputeIndexFromShapeAndSubsUsingSequence(shape,
                typename SequenceGenerator<sizeof...(Ds)>::type(), subs);
        }

        // ind2sub
        template <int ...Ds, class ...Subs>
        inline void SubsFromIndex(const Shape<Ds...> & shape, int index, Subs &... subs){
            static_assert(sizeof...(Ds) == sizeof...(Subs), "subcript length not compatible with shape rank!");
            ComputeSubsFromShapeAndIndexUsingSequence(shape, index,
                typename SequenceGenerator<sizeof...(Ds)>::type(), std::tie(subs...));
        }

        template <int ...Ds>
        inline void SubsFromIndex(const Shape<Ds...> & shape, int index, std::array<int, sizeof...(Ds)> & subs){
            ComputeSubsFromShapeAndIndexUsingSequence(shape, index,
                typename SequenceGenerator<sizeof...(Ds)>::type(), subs);
        }


        namespace {

            ////////////////////////
            template <int ...ADs, int ...BDs, int ...S>
            inline Shape<ADs...> MakeShapeFromDynamicSizesUsingSequence(const Shape<BDs...> & from, Sequence<S...>) {
                return Shape<ADs...>(from.dsize<S>()...);
            }

            template <int ...ADs, int ...BDs, int ...AS, int ...BS>
            inline Shape<ADs..., BDs...> MakeDescartesProductShapeUsingSequence(
                const Shape<ADs...> & froma, const Shape<BDs...> & fromb,
                Sequence<AS...>, Sequence<BS...>) {
                return Shape<ADs..., BDs...>(froma.dsize<AS>()..., fromb.dsize<BS>()...);
            }

            template <int M, int N, int P, int ...ADs, int ...BDs, int ...AS, int ...BS>
            inline Shape<M, P, ADs..., BDs...> MakeMatrixProductShapeUsingSequence(
                const Shape<M, N, ADs...> & a, const Shape<N, P, BDs...> & b,
                Sequence<AS...>, Sequence<BS...>) {
                return Shape<M, P, ADs..., BDs...>(a.dsize<AS>()..., b.dsize<BS>()...);
            }

        }

        // make a new shape from dynamic sizes of an existing shape
        template <int ...ADs, int ...BDs>
        inline Shape<ADs...> MakeShapeFromDynamicSizes(const Shape<BDs...> & from) {
            using FromShapeType = Shape<BDs...>;
            using ToShapeType = Shape<ADs...>;
            static_assert(FromShapeType::DynamicNum == ToShapeType::DynamicNum, "dynamic sizes count not compitable!");
            return MakeShapeFromDynamicSizesUsingSequence<ADs...>(from,
                typename SequenceGenerator<FromShapeType::DynamicNum>::type());
        }

        // make a new shape to represent descartes product
        template <int ...ADs, int ...BDs>
        inline Shape<ADs..., BDs...> MakeDescartesProductShape(const Shape<ADs...>& a, const Shape<BDs...>& b) {
            static const int DNA = Shape<ADs...>::DynamicNum;
            static const int DNB = Shape<BDs...>::DynamicNum;
            return MakeDescartesProductShapeUsingSequence(a, b,
                typename SequenceGenerator<DNA>::type(),
                typename SequenceGenerator<DNB>::type());
        }





        // dynamic shape
        class DShape {
        public:
            inline DShape() {}
            template <int ... Ds>
            inline DShape(const Shape<Ds...> & s) : _lengths(MakeVectorFromShape(s)) {}

            inline int rank() const { return _lengths.size(); }
            inline int size(int ord) const { return _lengths[ord]; }
            inline int volume() const { int v = 1; for (int l : _lengths){ v *= l; } return v; }
            inline void resize(int ord, int sz) { _lengths[ord] = sz; }

        private:
            std::vector<int> _lengths;
        };

        */





        class ExpressionGraph;
        struct Op;

        using GraphType = Mesh<std::shared_ptr<Op>>;
        using EHandle = GraphType::VertHandle;
        using CHandle = GraphType::HalfHandle;

        template <class DataT>
        struct OpWithAValue;
        
        // operator base
        struct Op {

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

            virtual std::shared_ptr<Op> makeZero() const = 0;
            virtual std::shared_ptr<Op> makeOne() const = 0;
            virtual std::shared_ptr<Op> makePlus() const = 0; // support > 2 args
            virtual std::shared_ptr<Op> makeMult() const = 0; // support > 2 args

            template <class DataT>
            const OpWithAValue<DataT> & as() const { return *static_cast<const OpWithAValue<DataT> *>(this); }

        };



        template <class DataT>
        struct OpWithAValue : public Op {
            virtual DataT value() const = 0;
            virtual std::shared_ptr<Op> makeZero() const override { return std::make_shared<OpWithAZeroValue<DataT>>(); }
            virtual std::shared_ptr<Op> makeOne() const override { return std::make_shared<OpWithAOneValue<DataT>>(); }
            virtual std::shared_ptr<Op> makePlus() const override { return std::make_shared<Plus<DataT>>(); }
            virtual std::shared_ptr<Op> makeMult() const override { return std::make_shared<Mult<DataT>>(); }
        };

        template <class DataT>
        struct OpWithAZeroValue : public OpWithAValue<DataT> {
            virtual DataT value() const override{ return DataTraits<DataT>::zero(); }
        };

        template <class DataT>
        struct OpWithAOneValue : public OpWithAValue<DataT> {
            virtual DataT value() const override{ return DataTraits<DataT>::one(); }
        };

        template <class DataT>
        struct OpWithACache : public OpWithAValue<DataT> {
            inline OpWithACache(const DataT & d = DataTraits<DataT>::zero()) : cache(d){}
            virtual DataT value() const override { return cache; }
            
            virtual void eval(ExpressionGraph const * const g, std::vector<EHandle> && inputs) override {}
            virtual std::vector<EHandle> makeInputDerivatives(ExpressionGraph * const g,
                EHandle self,
                std::vector<EHandle> && inputs,
                EHandle outputDerivsSum) const override {
                return std::vector<EHandle>();
            }

            DataT cache;
        };



        template <class DataT>
        struct Plus : public OpWithACache<DataT> {
            // evaluate and get result
            virtual void eval(ExpressionGraph const * const g, std::vector<EHandle> && inputs) {
                cache = DataTraits<DataT>::zero();
                for (auto h : inputs){
                    cache += g->op(h)->as<DataT>().value();
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


        template <class DataT>
        struct Mult : public OpWithACache<DataT> {
            // evaluate and get result
            virtual void eval(ExpressionGraph const * const g, std::vector<EHandle> && inputs) {
                cache = DataTraits<DataT>::one();
                for (auto h : inputs){
                    cache *= g->op(h)->as<DataT>().value();
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
                        g->addMult(toMults);
                }

                return inputDerivs;
            }
        };




        class ExpressionGraph {
        public:
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

            // get expression op
            inline const Op * op(EHandle h) const { return _g.data(h).get(); }
            template <class DataT>
            inline DataT value(EHandle h) const { return op(h)->as<DataT>().value(); }

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
            EHandle addPlus(const std::vector<EHandle>& inputs) { return addNode(op(inputs.front())->makePlus(), inputs); }
            EHandle addMult(const std::vector<EHandle>& inputs) { return addNode(op(inputs.front())->makeMult(), inputs); }
            template <class DataT>
            EHandle addValue(const DataT & d) { return addNode(std::make_shared<OpWithACache<DataT>>(d)); }

            // evaluate the expression
            void evaluate(EHandle h) const;
            template <class DataT>
            inline DataT evaluate(EHandle h) const { evaluate(h);  return op(h)->as<DataT>().value(); }

            // create derivative graph
            // returns the derivative nodes of vars
            std::vector<EHandle> createDerivatives(EHandle cost, const std::vector<EHandle> & vars);
            template <class ... EHandleTs>
            inline std::tuple<EHandleTs...> createDerivatives(EHandle cost, EHandleTs... vars){
                return createDerivativesUsingSequence(cost, SequenceGenerator<sizeof...(EHandleTs)>::type(), vars...);
            }

        private:
            template <class ... EHandleTs, int ... S>
            inline std::tuple<EHandleTs...> createDerivativesUsingSequence(EHandle cost, Sequence<S...>, EHandleTs... vars){
                std::vector<EHandle> derivs = createDerivatives(cost, std::vector<EHandle>{vars...});
                return std::make_tuple(derivs[S]...);
            }

            //// OpT must be a subclass of Op
            //// typename OpT::ElementType MUST exists
            //template <class OpT, int ...Ds, class ...OpArgTs>
            //inline Expr<typename OpT::ElementType, Ds...> addExpression(const Shape<Ds...> & shape, 
            //    std::initializer_list<EHandle> inputs, OpArgTs... opArgs){
            //    return Expr<typename OpT::ElementType, Ds...>(this, addNode(shape, std::make_shared<OpT>(opArgs...), inputs), shape);
            //}
           

        private:
            GraphType _g;
        };

        


        


       



        /*
        template <class ET, int ...Ds>
        struct Expr {
            inline Expr(typename ExpressionGraph<ET> * const g, EHandle<ET> h, const Shape<Ds...> & s)
                : graph(g), handle(h), resultShape(s){}

            template <class ...Ints>
            inline ET operator()(Ints... ints) const { 
                return graph->evaluate(handle, IndexFromSubs(resultShape, ints...)); 
            }

            Shape<Ds...> resultShape;
            ExpressionGraph<ET> const * const graph;
            EHandle<ET> handle;
        };


        template <class ET>
        struct Cache {
            inline Cache(int len) : data(new ET[len]), needsUpdate(true){}
            inline ~Cache(){ delete[] data; }
            ET * data; 
            bool needsUpdate;
        };




        // operator without inputs
        template <class ET>
        struct SrcOp : public Op<ET> {
            inline SrcOp(bool isConst = false) : Op<ET>(isConst){}

            virtual std::vector<EHandle<ET>> makeInputDerivatives(ExpressionGraph<ET> * const g,
                EHandle<ET> self,
                std::vector<EHandle<ET>> && inputs,
                EHandle<ET> outputDerivsSum) const {

                assert(inputs.size() == 0);
                return std::vector<EHandle<ET>>();
            }
        };

        // constant
        template <class ET>
        struct Const : public SrcOp<ET> {
            inline Const(const ET & c) : SrcOp<ET>(true), constant(c){}
            
            virtual ET eval(ExpressionGraph<ET> const * const g, 
                std::vector<EHandle<ET>> && inputs,
                EHandle<ET> self, 
                int index) const {

                return constant;
            }

            const ET constant;
        };

        // variable
        template <class ET>
        struct Var : public SrcOp<ET> {
            inline Var(const ET * c) : SrcOp<ET>(false), value(c){}

            virtual ET eval(ExpressionGraph<ET> const * const g, 
                std::vector<EHandle<ET>> && inputs,
                EHandle<ET> self, int index) const {
                return value[index];
            }

            const ET * value;
        };

        // elementwise plus
        template <class ET>
        struct Plus : public Op<ET> {
            virtual ET eval(ExpressionGraph<ET> const * const g,
                std::vector<EHandle<ET>> && inputs,
                EHandle<ET> self, int index) const {
                ET s = 0;
                for (auto i : inputs){
                    s += g->evaluate(i, index);
                }
                return s;
            }

            virtual std::vector<EHandle<ET>> makeInputDerivatives(ExpressionGraph<ET> * const g,
                EHandle<ET> self,
                std::vector<EHandle<ET>> && inputs,
                EHandle<ET> outputDerivsSum) const {

                std::vector<EHandle<ET>> inputDerivs(inputs.size(), outputDerivsSum);
                return inputDerivs;
            }            
        };

        // elementwise mul
        template <class ET>
        struct Mul : public Op<ET> {            
            virtual ET eval(ExpressionGraph<ET> const * const g,
                std::vector<EHandle<ET>> && inputs,
                EHandle<ET> self,
                int index) const {

                ET s = 1;
                for (auto i : inputs){
                    s *= g->evaluate(i, index);
                }
                return s;
            }

            virtual std::vector<EHandle<ET>> makeInputDerivatives(ExpressionGraph<ET> * const g,
                EHandle<ET> self,
                std::vector<EHandle<ET>> && inputs,
                EHandle<ET> outputDerivsSum) const {

                std::vector<EHandle<ET>> inputDerivs(inputs.size());
                for (int i = 0; i < inputs.size(); i++){
                    std::vector<EHandle<ET>> toMults(1, outputDerivsSum);
                    toMults.reserve(inputs.size());
                    for (int j = 0; j < inputs.size(); j++){
                        if (j == i)
                            continue;
                        toMults.push_back(inputs[j]);
                    }
                    inputDerivs[i] = toMults.size() == 1 ? toMults.front() :
                        g->addNode<Mul<ET>>(g->data(self).dshape, toMults);
                }

                return inputDerivs;
            }
        };

        // elementwise exponent
        template <class ET>
        struct Exp : public Op<ET> {
            
            virtual ET eval(ExpressionGraph<ET> const * const g,
                std::vector<EHandle<ET>> && inputs,
                EHandle<ET> self,
                int index) const {

                assert(inputs.size() == 1);
                return exp(g->evaluate(inputs.front(), index));                
            }

            virtual std::vector<EHandle<ET>> makeInputDerivatives(ExpressionGraph<ET> * const g,
                EHandle<ET> self,
                std::vector<EHandle<ET>> && inputs,
                EHandle<ET> outputDerivsSum) const {
           
                return{ g->addNode<Mul<ET>>(g->data(self).dshape,
                    outputDerivsSum,
                    g->addNode<Exp<ET>>(g->data(self).dshape)) };
            }

        };

        // matrix mul
        template <class ET>
        struct MatrixMul : public Op<ET> {

            virtual ET eval(ExpressionGraph<ET> const * const g,
                std::vector<EHandle<ET>> && inputs,
                EHandle<ET> self,
                int index) const {

                assert(inputs.size() == 2);
                
                // TODO
                return exp(g->evaluate(inputs.front(), index));
            }

            virtual std::vector<EHandle<ET>> makeInputDerivatives(ExpressionGraph<ET> * const g,
                EHandle<ET> self,
                std::vector<EHandle<ET>> && inputs,
                EHandle<ET> outputDerivsSum) const {

                // TODO
                return{ g->addNode<Mul<ET>>(g->data(self).dshape,
                    outputDerivsSum,
                    g->addNode<Exp<ET>>(g->data(self).dshape)) };
            }

        };







        // element wise functor
        //namespace fun {

        //    /// zeroary functor

        //    // const rational number
        //    template <class ET, int N, int D = 1>
        //    struct RationalConst {
        //        inline ET operator()() const { return N / D; }
        //    };

        //    template <class ET>
        //    using Zero = RationalConst<ET, 0>;
        //    template <class ET>
        //    using One = RationalConst<ET, 1>;

        //    template <class ET, int N, int D>
        //    Zero<ET> MakeGradientFunctor(RationalConst<ET, N, D>) { return Zero<ET>(); }

        //    // const number
        //    template <class ET>
        //    struct Constant {
        //        inline Constant(const ET & c) : constant(c){}
        //        inline ET operator()() const { return constant; }
        //        ET constant;
        //    };

        //    template <class ET>
        //    Zero<ET> MakeGradientFunctor(Constant<ET>){ return Zero<ET>(); }

        //    
        //    /// unary functor

        //    // + scalar
        //    template <class ET>
        //    struct PlusScalar {
        //        inline PlusScalar(const ET & s) : scalar(s){}
        //        inline ET operator()(const ET & a) const {
        //            return a + scalar;
        //        }
        //        ET scalar;
        //    };

        //    template <class ET>
        //    One<ET> MakeGradientFunctor(PlusScalar<ET> fun){
        //        return One<ET>();
        //    }

        //    // * scalar
        //    template <class ET>
        //    struct MultScalar {
        //        inline MultScalar(const ET & s) : scalar(s){}
        //        inline ET operator()(const ET & a) const {
        //            return a * scalar;
        //        }
        //        ET scalar;
        //    };

        //    template <class ET>
        //    Constant<ET> MakeGradientFunctor(MultScalar<ET> fun){
        //        return Constant<ET>(fun.scalar);
        //    }

        //    // exp
        //    template <class ET>
        //    struct Exp {
        //        inline ET operator()(const ET & a) const {
        //            return exp(a);
        //        }
        //    };

        //    template <class ET>
        //    Exp<ET> MakeGradientFunctor(Exp<ET> fun){
        //        return fun;
        //    }

        //}


 





        /*


        template <class ET, int ...A, int ...B>
        Expr<ET, Shape<A...>::StaticVolume, Shape<B...>::StaticVolume> JacobianMatrix(const Expr<ET, A...> & a, const Expr<ET, B...> & b) {

            throw "";
        }

        template <class ET, int ...A, int ...B>
        Expr<ET, Shape<A...>::StaticVolume, Shape<B...>::StaticVolume, Shape<B...>::StaticVolume> HessianMatrix(const Expr<ET, A...> & a, const Expr<ET, B...> & b) {

            throw "";
        }

        template <class ET, int A, int B, int C>
        Expr<ET, A, C> MatrixMult(const Expr<ET, A, B> & a, const Expr<ET, B, C> & b) {

            throw "";
        }


        template <class ET, int ... Ds>
        Expr<ET, Ds...> operator - (const Expr<ET, Ds...> & e) {
            throw "";
        }

        template <class ET, int ... Ds>
        Expr<ET, Ds...> exp(const Expr<ET, Ds...> & e) {
            throw "";
        }

        template <class ET, int ... Ds>
        Expr<ET, Ds...> operator + (const Expr<ET, Ds...> & e, const ET & scalar) {
            throw "";
        }

        template <class ET, int ... Ds>
        Expr<ET, Ds...> operator + (const ET & scalar, const Expr<ET, Ds...> & e) {
            throw "";
        }

        template <class ET, int ... Ds>
        Expr<ET, Ds...> operator + (const Expr<ET, Ds...> & e1, const Expr<ET, Ds...> & e2) {
            throw "";
        }

        template <class ET, int ... Ds>
        Expr<ET, Ds...> operator - (const Expr<ET, Ds...> & e, const ET & scalar) {
            throw "";
        }

        template <class ET, int ... Ds>
        Expr<ET, Ds...> operator - (const ET & scalar, const Expr<ET, Ds...> & e) {
            throw "";
        }

        template <class ET, int ... Ds>
        Expr<ET, Ds...> operator - (const Expr<ET, Ds...> & e1, const Expr<ET, Ds...> & e2) {
            throw "";
        }

        template <class ET, int ... Ds>
        Expr<ET, Ds...> operator * (const Expr<ET, Ds...> & e, const ET & scalar) {
            throw "";
        }

        template <class ET, int ... Ds>
        Expr<ET, Ds...> operator * (const ET & scalar, const Expr<ET, Ds...> & e) {
            throw "";
        }

        template <class ET, int ... Ds>
        Expr<ET, Ds...> operator * (const Expr<ET, Ds...> & e1, const Expr<ET, Ds...> & e2) {
            throw "";
        }

        template <class ET, int ... Ds>
        Expr<ET, Ds...> operator / (const Expr<ET, Ds...> & e, const ET & scalar) {
            throw "";
        }

        template <class ET, int ... Ds>
        Expr<ET, Ds...> operator / (const ET & scalar, const Expr<ET, Ds...> & e) {
            throw "";
        }

        template <class ET, int ... Ds>
        Expr<ET, Ds...> operator / (const Expr<ET, Ds...> & e1, const Expr<ET, Ds...> & e2) {
            throw "";
        }

        template <class ET, int ... Ds>
        Expr<ET> SumAll(const Expr<ET, Ds...> & e) {
            throw "";
        }




        template <class ET>
        struct Performer {

        };*/


    }
}




 
#endif