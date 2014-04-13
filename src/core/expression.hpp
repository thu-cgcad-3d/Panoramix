#ifndef PANORAMIX_CORE_EXPRESSION_HPP
#define PANORAMIX_CORE_EXPRESSION_HPP

#include <iostream>
#include <initializer_list>

#include "basic_types.hpp"
#include "template_utilities.hpp"

namespace panoramix {
    namespace core {

        using Subscript = std::initializer_list<int>;
        using SubscriptIter = std::initializer_list<int>::const_iterator;

        // tensor traits
        template <class T>
        struct TensorTraits {
            static const int Rank = 0;
            using ElementType = T;
            static ElementType ValueAt(const T & t, SubscriptIter subs) { return t; }
        };

        template <class T, class AllocT>
        struct TensorTraits<std::vector<T, AllocT>> {
            static const int Rank = TensorTraits<T>::Rank + 1;
            using ElementType = typename TensorTraits<T>::ElementType;
            static ElementType ValueAt(const std::vector<T, AllocT> & m, SubscriptIter subs) {
                return TensorTraits<T>::ValueAt(m[*subs], subs + 1); 
            }
        };

        template <class T, int M>
        struct TensorTraits<Vec<T, M>> {
            static const int Rank = TensorTraits<T>::Rank + 1;
            using ElementType = typename TensorTraits<T>::ElementType;
            static ElementType ValueAt(const Vec<T, M> & m, SubscriptIter subs) {
                return TensorTraits<T>::ValueAt(m(*subs), subs + 1);
            }
        };

        template <class T, int M>
        struct TensorTraits<T[M]> {
            static const int Rank = TensorTraits<T>::Rank + 1;
            using ElementType = typename TensorTraits<T>::ElementType;
            static ElementType ValueAt(const T(& m)[M], SubscriptIter subs) {
                return TensorTraits<T>::ValueAt(m[*subs], subs + 1);
            }
        };

        template <class T>
        struct TensorTraits<T[]> {
            static const int Rank = TensorTraits<T>::Rank + 1;
            using ElementType = typename TensorTraits<T>::ElementType;
            static ElementType ValueAt(const T(&m)[], SubscriptIter subs) {
                return TensorTraits<T>::ValueAt(m[*subs], subs + 1);
            }
        };

        template <class T, int M, int N>
        struct TensorTraits<Mat<T, M, N>> {
            static const int Rank = TensorTraits<T>::Rank + 2;
            using ElementType = typename TensorTraits<T>::ElementType;
            static ElementType ValueAt(const Mat<T, M, N> & m, SubscriptIter subs) {
                return TensorTraits<T>::ValueAt(m(*(subs), *(subs + 1)), subs + 1);
            }
        };


        
        // expression content structs

        // cache for expression calculation


        // base content
        template <class ElementT>
        struct ExpressionContent {
            using ElementType = ElementT;
            virtual ElementType eval(Subscript subs) const = 0;
            inline ElementType operator()(Subscript subs) const { return eval(subs); }
            virtual ElementType gradient(Subscript subs) const { return 0; }
        };

        template <class ElementT>
        using ExpressionContentPtr = std::shared_ptr<ExpressionContent<ElementT>>;

        namespace  { // details
            template <class TensorOperatorT, class TupleT, int ...S>
            inline auto TensorInvokeWithEachTupleArg(TensorOperatorT fun,
                TupleT args, Subscript subs, Sequence<S...>)
                -> decltype(fun(subs, *std::get<S>(args) ...)) {
                return fun(subs, *std::get<S>(args) ...);
            }
        }

        // tensor operation content
        template <class ElementT, class TensorOperatorT, class ...ArgTs>
        struct TensorExpressionContent : public ExpressionContent<ElementT> {
            inline TensorExpressionContent(TensorOperatorT f, ExpressionContentPtr<ArgTs>... a)
            : fun(f), args(std::forward_as_tuple(a...)) {}

            virtual ElementType eval(Subscript subs) const override {
                return TensorInvokeWithEachTupleArg(fun, args, subs,
                    typename SequenceGenerator<sizeof...(ArgTs)>::type());
            }

            std::tuple<ExpressionContentPtr<ArgTs>...> args;
            TensorOperatorT fun;
        };




        // expression class
        template <class ElementT>
        class Expression {
        public:
            Expression() :_p(nullptr){}
            explicit Expression(ExpressionContentPtr<ElementT> p)
                : _p(p) {}

            inline ExpressionContentPtr<ElementT> ptr() const { return _p; }
            inline ElementT eval(Subscript subs) const { return _p->eval(subs); }
            inline ElementT operator()() const { return eval({}); }
            inline ElementT operator()(int i) const { return eval({ i }); }
            inline ElementT operator()(int r, int c) const { return eval({ r, c }); }
            inline ElementT operator()(int r, int c, int d) const { return eval({ r, c, d }); }

            template <class DerivElementT>
            inline ElementT derivative(Subscript subs, 
                const Expression<DerivElementT> & d, Subscript dsubs) const;

        private:
            ExpressionContentPtr<ElementT> _p;
        };






        // operators
        template <class ResultT>
        struct TensorOperatorBase {
            using ResultType = ResultT;
        };


        // Constant
        template <class T>
        struct ConstantOperator : public TensorOperatorBase<typename TensorTraits<T>::ElementType> {
            inline ConstantOperator(){}
            inline ConstantOperator(const T & v) : value(v){}
            inline ConstantOperator(T && v) : value(v){}
            inline ResultType operator()(Subscript subs) const {
                return TensorTraits<T>::ValueAt(value, subs.begin());
            }
            T value;
        };

        template <class T>
        inline Expression<typename TensorTraits<T>::ElementType>
            MakeConstantValue(const T & value){
            using ResultType = typename TensorTraits<T>::ElementType;
            return Expression<ResultType>
                (std::make_shared<TensorExpressionContent<ResultType, ConstantOperator<T>>>
                (ConstantOperator<T>(value)));
        }

        template <class T>
        inline Expression<typename TensorTraits<std::remove_reference_t<T>>::ElementType>
            MakeConstantValue(T && value){
            using PureT = std::remove_reference_t<T>;
            using ResultType = typename TensorTraits<PureT>::ElementType;
            return Expression<ResultType>
                (std::make_shared<TensorExpressionContent<ResultType, ConstantOperator<PureT>>>
                (ConstantOperator<PureT>(value)));
        }


        // Variable
        template <class T>
        struct VariableRefOperator : public TensorOperatorBase<typename TensorTraits<T>::ElementType> {
            inline VariableRefOperator(T * v) : value(v){}
            inline ResultType operator()(Subscript subs) const {
                return TensorTraits<T>::ValueAt(*value, subs.begin());
            }
            T * value;
        };

        template <class T>
        inline Expression<typename TensorTraits<T>::ElementType> MakeVariableAt(T * v) {
            using ResultType = typename TensorTraits<T>::ElementType;
            return Expression<ResultType>
                (std::make_shared<TensorExpressionContent<ResultType, VariableRefOperator<T>>>
                (VariableRefOperator<T>(v)));
        }



        // Element Wise Operation
        template <class ResultT, class ElementWiseFunctorT, class ...ArgTs>
        struct ElementWiseOperator : public TensorOperatorBase<ResultT> {
            inline ElementWiseOperator(ElementWiseFunctorT ef) : eleFun(eleFun){}
            inline ResultType operator()(Subscript subs, const ExpressionContent<ArgTs> &... args) const {
                return eleFun(args.eval(subs)...);
            }
            
            ElementWiseFunctorT eleFun;
        };

        // perform element wise operation
        template <class ElementWiseFunctorT, class ...ArgTs>
        inline auto PerformElementWiseOperation(ElementWiseFunctorT fun, const Expression<ArgTs> &... args)
            ->Expression<decltype(fun(args()...))> {
            using ResultType = decltype(fun(args()...));
            using FunctorType = ElementWiseOperator<ResultType, ElementWiseFunctorT, ArgTs...>;
            return Expression<ResultType>
                (std::make_shared<TensorExpressionContent<ResultType, FunctorType, ArgTs...>>
                (FunctorType(fun), args.ptr()...));
        }

        

        

        







    }
}
 
#endif