#ifndef PANORAMIX_CORE_EXPRESSION_HPP
#define PANORAMIX_CORE_EXPRESSION_HPP

#include <iostream>
#include <utility>
#include <initializer_list>

#include "basic_types.hpp"
#include "template_utilities.hpp"
#include "mesh.hpp"

namespace panoramix {
    namespace core {

        enum DimensionType : int {
            Dynamic = 0
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
        class Shape<Dynamic, Ds...> : private Shape<Ds...>{
        public:
            using BaseType = Shape<Ds...>;
            static const bool LengthIsDynamic = true;
            static const int StaticLength = Dynamic;
            static const int StaticVolume = Dynamic;
            static const int Rank = sizeof...(Ds)+1;
            static const int DynamicNum = BaseType::DynamicNum + 1;

            // set size for Dynamic dimension
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

        // make a new shape to represent result of matrix multiplication
        template <int M, int N, int P, int ...ADs, int ...BDs>
        inline Shape<M, P, ADs..., BDs...> MakeMatrixProductShape(const Shape<M, N, ADs...> & a, const Shape<N, P, BDs...> & b) {
            static const int DNA = Shape<M, N, ADs...>::DynamicNum;
            static const int DNB = Shape<N, P, BDs...>::DynamicNum;
            return MakeMatrixProductShapeUsingSequence(a, b,
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



        template <class ET, int ...Ds>
        struct Expr;

        template <class ET>
        struct Op;

        template <class ET>
        struct Const;

        template <class ET>
        struct Plus;

        template <class ET>
        class ExpressionGraph {
        public:
            struct ExpressionData;
            struct ConnectionData;
            
            using GraphType = Mesh<ExpressionData, ConnectionData>;
            using ExpressionHandle = typename GraphType::VertHandle;
            using ConnectionHandle = typename GraphType::HalfHandle;
            
            struct ExpressionData {
                DShape dshape;
                std::shared_ptr<Op<ET>> op;
            };

            struct ConnectionData {
               
            };

        public:

            // checks whether this is a forward connection 
            // (from inputs to outputs, from old expressions to newly created expressions)
            inline bool isForwardConnection(ConnectionHandle h) const {
                auto & topo = _g.topo(h);
                return topo.from().id < topo.to().id;
            }
            inline bool isBackwardConnection(ConnectionHandle h) const { return !isForwardConnection(h); }

            // add expression
            // OpT must be a subclass of Op<ET>
            template <class OpT, class ...OpArgTs>
            ExpressionHandle addNode(const DShape & dshape, const std::vector<ExpressionHandle>& inputs, OpArgTs... opArgs) {
                auto h = _g.addVertex();
                _g.data(h).dshape = dshape;
                _g.data(h).op = std::make_shared<OpT>(opArgs...);
                for (auto & ih : inputs){
                    _g.addEdge(ih, h);
                }
                return h;
            }

            ET evaluate(ExpressionHandle h, int index) const {
                std::vector<ExpressionHandle> inputs;
                inputs.reserve(_g.topo(h).halfedges.size());
                for (auto hh : _g.topo(h).halfedges){
                    if (_g.exists(hh) && isBackwardConnection(hh)){
                        inputs.push_back(_g.topo(hh).to());
                    }
                }
                return _g.data(h).op->eval(this, std::move(inputs), h, index);
            }

            std::vector<ExpressionHandle> createDerivativeGraph(ExpressionHandle cost) {
                assert(_g.data(cost).dshape.rank() == 0); // cost MUST be a scalar!

                // the id of expression handles record the topological order
                // idToDeriv[exprHandle.id] stores the derivative handle of exprHandle
                std::vector<std::vector<ExpressionHandle>> idToDerivs(cost.id + 1);
                idToDerivs[cost.id] = { addNode<Const<ET>>(_g.data(cost).dshape, {}, 1) };

                std::vector<ExpressionHandle> idToDerivsSumTable(cost.id + 1);
                for (int i = cost.id; i >= 0; i--){
                    ExpressionHandle self(i);
                    if (!_g.exists(self))
                        continue;

                    // this expression is not used by cost at all
                    if (idToDerivs[self.id].empty())
                        continue;

                    // sum all derivatives of self outputs
                    auto idToDerivsSum = idToDerivs[self.id].size() == 1 ? idToDerivs[self.id].front() :
                        addNode<Plus<ET>>(_g.data(self).dshape, idToDerivs[self.id]);

                    // get all inputs
                    std::vector<ExpressionHandle> inputs;
                    for (auto conh : _g.topo(self).halfedges){
                        if (_g.exists(conh) && isBackwardConnection(conh)){
                            auto input = _g.topo(conh).to();
                            inputs.push_back(input);
                        }                        
                    }

                    // compute input derivatives
                    std::vector<ExpressionHandle> inputDerivs =
                        _g.data(self).op->makeInputDerivatives(this, self,
                        std::move(inputs), idToDerivsSum);

                    assert(inputDerivs.size() == inputs);

                    // store the input derivatives
                    for (int k = 0; k < inputs.size(); k++){
                        idToDerivs[inputs[k].id].push_back(inputDerivs[k]);
                    }

                    idToDerivsSumTable[self.id] = idToDerivsSum;
                } 

                return idToDerivsSumTable;
            }



        public:
            template <class OpT, int ...Ds, class ...OpArgTs>
            inline Expr<ET, Ds...> addExpression(const Shape<Ds...> & shape, 
                std::initializer_list<ExpressionHandle> inputHandles, OpArgTs... opArgs){
                return Expr<ET, Ds...>(this, addNode<OpT, OpArgTs...>(shape, inputHandles, opArgs...), shape);
            }

            inline const GraphType & graph() const { return _g; }

        private:
            GraphType _g;
        };


        template <class ET, int ...Ds>
        struct Expr {
            inline Expr(typename ExpressionGraph<ET> * const g,
                typename ExpressionGraph<ET>::ExpressionHandle h,
                const Shape<Ds...> & s)
                : graph(g), handle(h), resultShape(s){}

            template <class ...Ints>
            inline ET eval(Ints... ints) const { 
                return graph->evaluate(handle, IndexFromSubs(resultShape, ints...)); 
            }

            Shape<Ds...> resultShape;
            typename ExpressionGraph<ET> const * const graph;
            typename ExpressionGraph<ET>::ExpressionHandle handle;
        };


        // operator base
        template <class ET>
        struct Op {
            using ExpressionData = typename ExpressionGraph<ET>::ExpressionData;
            using ExpressionHandle = typename ExpressionGraph<ET>::ExpressionHandle;

            inline Op(bool isConst = false) : isConstant(isConst){}

            // evaluate result
            virtual ET eval(ExpressionGraph<ET> const * const g, 
                std::vector<ExpressionHandle> && inputs,
                ExpressionHandle self, int index) const = 0;
            
            // make input derivative expressions based on inputs and output derivative expressions 
            virtual std::vector<ExpressionHandle> makeInputDerivatives(ExpressionGraph<ET> * const g, 
                ExpressionHandle self,
                std::vector<ExpressionHandle> && inputs,
                ExpressionHandle outputDerivsSum) const {
                return std::vector<ExpressionHandle>(inputs.size());
            }

            const bool isConstant;
        };

        // constant
        template <class ET>
        struct Const : public Op<ET> {
            inline Const(const ET & c) : Op<ET>(true), constant(c){}
            
            virtual ET eval(ExpressionGraph<ET> const * const g, 
                std::vector<ExpressionHandle> && inputs,
                ExpressionHandle self, int index) const {
                return constant;
            }

            const ET constant;
        };

        // variable
        template <class ET>
        struct Var : public Op<ET> {
            inline Var(const ET * c) : Op<ET>(), value(c){}

            virtual ET eval(ExpressionGraph<ET> const * const g, 
                std::vector<ExpressionHandle> && inputs,
                ExpressionHandle self, int index) const {
                return value[index];
            }

            const ET * value;
        };


        template <class ET>
        struct Plus : public Op<ET> {

            virtual ET eval(ExpressionGraph<ET> const * const g,
                std::vector<ExpressionHandle> && inputs,
                ExpressionHandle self, int index) const {
                ET s = 0;
                for (auto i : inputs){
                    s += g->evaluate(i, index);
                }
                return s;
            }

            virtual std::vector<ExpressionHandle> makeInputDerivatives(ExpressionGraph<ET> * const g,
                ExpressionHandle self,
                std::vector<ExpressionHandle> && inputs,
                ExpressionHandle outputDerivsSum) const {

                return std::vector<ExpressionHandle>(inputs.size());
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

        //template <class ET>
        //struct Var;
        //template <class ET>
        //struct Node;


        //template <class ET>
        //struct ExprGraph {
        //    struct Vertex {
        //        bool isVar;
        //        Node<ET> * node;
        //        Var<ET> * var;
        //    };

        //    using Graph = Mesh<Vertex>;
        //    Graph g;
        //};

        //template <class ET>
        //using VarHandle = typename ExprGraph<ET>::Graph::HalfHandle;
        //using NodeHandle = typename ExprGraph<ET>::Graph::VertHandle;
        //

        //
        //// vars
        //template <class ET>
        //struct Var {
        //    virtual const ET * valueArray(const ExprGraph<ET> & g, VarHandle<ET> h) const { return nullptr; }
        //    virtual ET valueAt(int index, const ExprGraph<ET> & g, VarHandle<ET> h) const { return 0; }
        //    virtual int valueNum() const { return std::numeric_limit<int>::max(); }
        //};

        //template <class ET>
        //struct VarFinite : public Var<ET> {
        //    inline VarFinite(int n) : num(n){}
        //    virtual int valueNum() const { return num; }
        //    int num;
        //};

        //template <class ET, class ShapeT>
        //struct VarShaped : public VarFinite<ET> {
        //    inline VarShaped(const ShapeT & s) : VarFinite<ET>(s.volume()), shape(s){}
        //    ShapeT shape;
        //};

        //template <class ET, class ShapeT>
        //struct Constant : public VarShaped<ET, ShapeT> {
        //    inline Constant(const ET & c, const ShapeT & s) : VarShaped<ET, ShapeT>(s), constant(c){}
        //    virtual ET valueAt(int index, const ExprGraph<ET> & g, VarHandle<ET> h) const { return constant; }
        //    ET constant;
        //};
        //
        //// nodes
        //template <class ET>
        //struct Node {

        //    // back propagate
        //    //virtual void createInputGradients(
        //    //    const std::vector<Var<ET>*> & outputGrads, 
        //    //     std::vector<Var<ET>*> & inputGrads) {}
        //    
        //    virtual 
        //    virtual Node<ET>* createBackPropagateNode() const { return nullptr; }
        //};




        //





        //// convenient base for operators
        //template <class ET, class ShapeOutT, class ...ShapeInTs>
        //struct OpBase {
        //    using ShapeOutType = ShapeOutT;

        //    std::tuple<Var<ET, ShapeInTs> *...> inputs;
        //    Var<ET, ShapeOutT> * output;
        //};

        //template <class ET, class ShapeT>
        //struct PlusOp : public OpBase<ET, ShapeT, ShapeT, ShapeT>{

        //    inline void grad()

        //};


        //template <class ET>
        //class ExprCache;
        //
        //// expression tools
        //template <class ET>
        //struct ExprPerformer {
        //    virtual int shapeVolume() const { return 1; }
        //    virtual ET value(int index) const { return 0; }
        //    virtual ET value(int index, ExprCache<ET> * cache) const { return value(index); }
        //    virtual ExprPerformer * createGradientPerformer(const ExprPerformer * d) const { return nullptr; }
        //};

        //// square eye performer
        //template <class ET>
        //struct SquareEyeConstPerformer : public ExprPerformer<ET> {
        //    inline SquareEyeConstPerformer(int sv1, int sv2) : inShapeVolume1(sv1), inShapeVolume2(sv2){}
        //    virtual int shapeVolume() const { return inShapeVolume1 * inShapeVolume2; }
        //    virtual ET value(int index) const {
        //        return (index % inShapeVolume1 == index / inShapeVolume1) ? 1 : 0;
        //    }
        //    virtual ExprPerformer * createGradientPerformer(const ExprPerformer * d) const {
        //        return this == d ? new SquareEyeConstPerformer(shapeVolume(), d->shapeVolume()) 
        //            : new ZeroPerformer<ET>(shapeVolume() * d->shapeVolume());
        //    }
        //    int inShapeVolume1, inShapeVolume2;
        //};

        //// element wise plus performer
        //template <class ET>
        //struct ElementWisePlusPerformer : public ExprPerformer<ET> {
        //    inline ElementWisePlusPerformer(ExprPerformer * mm1, ExprPerformer * mm2)
        //        : m1(mm1), m2(mm2) {}
        //    virtual int shapeVolume() const { return std::min(m1->shapeVolume(), m2->shapeVolume()); }
        //    virtual ET value(int index) const { return m1->value(index) + m2->value(index); }
        //    virtual ExprPerformer * createGradientPerformer(const ExprPerformer * d) const {
        //        return this == d ? new SquareEyeConstPerformer<ET>(shapeVolume(), d->shapeVolume()) 
        //            : new ElementWisePlusPerformer<ET>(
        //            m1->createGradientPerformer(d), 
        //            m2->createGradientPerformer(d));
        //    }

        //    ExprPerformer * m1;
        //    ExprPerformer * m2;
        //};

        //// element wise mult performer
        //template <class ET>
        //struct ElementWiseMultPerformer : public ExprPerformer<ET> {
        //    inline ElementWiseMultPerformer(ExprPerformer * mm1, ExprPerformer * mm2) 
        //        : m1(mm1), m2(mm2){}
        //    virtual int shapeVolume() const { return std::min(m1->shapeVolume(), m2->shapeVolume()); }
        //    virtual ET value(int index) const { return m1->value(index) * m2->value(index); }
        //    virtual ExprPerformer * createGradientPerformer(const ExprPerformer * d) const {
        //        return this == d ? new SquareEyeConstPerformer<ET>(shapeVolume(), d->shapeVolume()) 
        //            : new ElementWisePlusPerformer<ET>(
        //            new ElementWiseMultPerformer(m1->createGradientPerformer(d), m2), 
        //            new ElementWiseMultPerformer(m1, m2->createGradientPerformer(d)));
        //    }

        //    ExprPerformer * m1;
        //    ExprPerformer * m2;
        //};


        //// zeroary performer [LEAF]
        //template <class ET, class FunT>
        //struct ElementWiseZeroaryPerformer : public ExprPerformer<ET> {
        //    inline ElementWiseZeroaryPerformer(int sv, FunT f = FunT()) : outShapeVolume(sv), fun(f){}
        //    virtual int shapeVolume() const { return outShapeVolume; }
        //    virtual ET value(int index) const { return fun(); }
        //    virtual ExprPerformer * createGradientPerformer(const ExprPerformer * d) const {
        //        return this == d ? new SquareEyeConstPerformer(shapeVolume(), d->shapeVolume())
        //            : new ZeroPerformer<ET>(shapeVolume() * d->shapeVolume());
        //    }

        //    FunT fun;
        //    int outShapeVolume;
        //};


        //// general unary performer on each single element
        //template <class ET, class FunT>
        //struct ElementWiseUnaryPerformer : public ExprPerformer<ET> {
        //    inline ElementWiseUnaryPerformer(ExprPerformer * mm, FunT f = FunT()) : m(mm), fun(f){}
        //    virtual int shapeVolume() const { return m->shapeVolume(); }
        //    virtual ET value(int index) const { return fun(m->value(index)); }
        //    virtual ExprPerformer * createGradientPerformer(const ExprPerformer * d) const {
        //        if (this == d)
        //            return new SquareEyeConstPerformer(shapeVolume(), d->shapeVolume());
        //        auto gradientFun = MakeGradientFunctor(fun);
        //        using GradientType = decltype(gradientFun);
        //        return new ElementWiseMultPerformer<ET>(
        //            new ElementWiseUnaryPerformer<GradientType>(m, gradientFun), 
        //            m->createGradientPerformer(d));
        //    }

        //    ExprPerformer * m;
        //    FunT fun;
        //};

        //// all sum
        ////template <class ET>
        ////struct AllSumPerformer : public ExprPerformer<ET> {
        ////    inline AllSumPerformer(ExprPerformer * mm) : m(mm){}
        ////    virtual int shapeVolume() const { return 1; }
        ////    virtual ET value(int index) const {
        ////        ET sum = 0;
        ////        for (int i = 0; i < shapeVolume; i++){
        ////            sum += m->value(i);
        ////        }
        ////        return sum;
        ////    }
        ////    virtual ExprPerformer * createGradientPerformer(const ExprPerformer * d) const {
        ////        if (this == d)
        ////            return new SquareEyeConstPerformer(shapeVolume(), d->shapeVolume());
        ////        return new OnePerformer<ET>(shapeVolume() * d->shapeVolume());
        ////    }

        ////    ExprPerformer * m;
        ////};


        //template <class ET, class Shape1, class Shape2>
        //struct MatrixMultPerformer : public ExprPerformer<ET> {
        //    static_assert(Shape1::Rank >= 2 && Shape2::Rank >= 2, "shape ranks must be over 2!");
        //    static const int Size10 = Shape1::StaticLength;
        //    static const int Size11 = Shape1::BaseType::StaticLength;
        //    static const int Size20 = Shape2::StaticLength;
        //    static const int Size21 = Shape2::BaseType::StaticLength;
        //    static_assert(Size11 == Size20, "shape sizes not compatible!");
        //    using ResultShapeType = decltype(MakeMatrixProductShape(std::declval<Shape1>(), std::declval<Shape2>()));

        //    MatrixMultPerformer(ExprPerformer * mm1, ExprPerformer * mm2, const Shape1 & s1, const Shape2 & s2)
        //        : m1(mm1), m2(mm2), shape1(s1), shape2(s2), resultShape(MakeMatrixProductShape(s1, s2)){
        //        assert(s1.size<1>() == s2.size<0>());
        //    }

        //    virtual int shapeVolume() const { return resultShape.volume(); }
        //    virtual ET value(int index) const {
        //        // convert index back to subscript
        //        std::array<int, ResultShapeType::Rank> resultSubs;
        //        SubsFromIndex(resultShape, index, resultSubs);
        //        // get subscripts of the two inputs
        //        std::array<int, Shape1::Rank> subs1;
        //        std::array<int, Shape2::Rank> subs2;
        //        std::copy(resultSubs.begin() + 3, resultSubs.end(), subs1.begin() + 3);
        //        std::copy(resultSubs.begin() + 3, resultSubs.end(), subs2.begin() + 3);
        //        subs1[0] = resultSubs[0];
        //        subs2[1] = resultSubs[1];
        //        static const int N = Size11;
        //        ET sum = 0;
        //        for (int i = 0; i < N; i++){
        //            subs1[1] = subs2[0] = i;
        //            sum += m1->value(IndexFromSubs(shape1, subs1)) * m2->value(IndexFromSubs(shape2, subs2));
        //        }
        //        return sum;
        //    }

        //    virtual ExprPerformer * createGradientPerformer(const ExprPerformer * d) const {
        //        
        //    }

        //    Shape1 shape1;
        //    Shape2 shape2;
        //    ResultShapeType resultShape;
        //    ExprPerformer * m1;
        //    ExprPerformer * m2;
        //};





        /*template <class ET, int ... Ds>
        class Expr {
        public:
            using ShapeType = Shape<Ds...>;

            template <class... Ints>
            explicit Expr(Ints ... ints) : _shape(ints...), _data(std::make_shared<ExprPerformer<ET>>()) {}

            Expr(std::shared_ptr<ExprPerformer<ET>> ptr, const ShapeType & s) : _shape(s), _data(ptr) {}

            inline ShapeType shape() const {
                return _shape;
            }

            Expr<ET, Ds..., 1> promote() const {
                return Expr<ET, Ds..., 1>(_data, MakeShapeFromDynamicSizes<Ds..., 1>(_shape));
            }

            template <int ... OtherDs>
            Expr<ET, OtherDs...> reshape(const Shape<OtherDs...> & s) const {
                static_assert(Shape<OtherDs...>::StaticVolume == ShapeType::StaticVolume, "shape volume not compatible!");
                assert(s.volume() == _shape.volume());
                return Expr<ET, OtherDs...>(_data, s);
            }

        private:
            ShapeType _shape;
            std::shared_ptr<ExprPerformer<ET>> _data;
        };


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