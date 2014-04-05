#ifndef PANORAMIX_CORE_EXPRESSION_HPP
#define PANORAMIX_CORE_EXPRESSION_HPP

#include <iostream>

#include "basic_types.hpp"
#include "template_utilities.hpp"

namespace panoramix {
    namespace core {

        namespace {

            template <int ... Dims>
            struct Dimension {};

            template <class T>
            struct DimensionOfType {
                using type = Dimension<1>;
            };

            template <class T, int M, int N>
            struct DimensionOfType<Mat<T, M, N>> {
                using type = Dimension<M, N>;
            };

            template <class A, class B>
            struct DescartesProductDimension {
                using type = Dimension<1>;
            };
            


            // fundamental structure
            // base of all expressions whose return type is T
            template <class T>
            struct ExpressionContentBase {
                virtual T eval() const = 0;
                virtual ExpressionContentBase * derivative(T * var) const { return nullptr; }
                virtual std::string texify() const { return ""; }
            };
            template <class T>
            using ExpressionContentBasePtr = std::shared_ptr<ExpressionContentBase<T>>;

            template <class FunctorT, class TupleT, int ...S>
            inline auto InvokeWithEachTupleArgEvaled(FunctorT fun, TupleT args, Sequence<S...>)
                -> decltype(fun(std::get<S>(args)->eval() ...)) {
                return fun(std::get<S>(args)->eval() ...);
            }

            template <class T, class FunctorT, class ...ArgTs>
            struct ExpressionContent : public ExpressionContentBase<T> {
                inline ExpressionContent(FunctorT f, ExpressionContentBasePtr<ArgTs>... a)
                : fun(f), args(std::forward_as_tuple(a...)) {}

                virtual T eval() const override {
                    return InvokeWithEachTupleArgEvaled(fun, args,
                        typename SequenceGenerator<sizeof...(ArgTs)>::type());
                }

                std::tuple<ExpressionContentBasePtr<ArgTs>...> args;
                FunctorT fun;
            };

            template <class T, class FunctorT, class ...ArgTs>
            inline ExpressionContentBasePtr<T> MakeExpressionPtr(FunctorT && fun,
                ExpressionContentBasePtr<ArgTs> && ... args) {
                return std::make_shared<ExpressionContent<T, FunctorT, ArgTs...>>(fun, args...);
            }

        }


        // operator structures

        // Constant
        template <class T>
        struct Constant {
            inline Constant(){}
            inline Constant(const T & v) : value(v){}
            inline T operator()(void) const {
                return value;
            }
            T value;
        };

        // Variable
        template <class T>
        struct VariableRef {
            inline VariableRef(T * v) : value(v){}
            inline T operator()(void) const {
                return *value;
            }
            T * value;
        };


        // the expression class
        template <class T>
        class Expression {
        public:
            Expression() :_p(nullptr){}
            Expression(const T & value)
                : _p(MakeExpressionPtr<T>(Constant<T>(value))) {}
            Expression(T * value)
                : _p(MakeExpressionPtr<T>(VariableRef<T>(value))) {}
            explicit Expression(ExpressionContentBasePtr<T> p)
                : _p(p) {}
            
            inline ExpressionContentBasePtr<T> ptr() const { return _p; }
            inline T eval() const { return _p->eval(); }

        private:
            ExpressionContentBasePtr<T> _p;
        };

        template <class T>
        inline Expression<T> MakeConstantValue(const T & value){ 
            return Expression<T>(value); 
        }
        
        template <class T>
        inline Expression<T> MakeVariableAt(T * v) {
            return Expression<T>(v);
        }


        // perform expression operation on existing functors
        template <class FunctorT, class ...ArgTs>
        inline auto PerformOperation(FunctorT fun, const Expression<ArgTs> &... args) 
            -> Expression<decltype(fun(args.eval()...))> {
            using ResultType = decltype(fun(args.eval()...));
            return Expression<ResultType>(MakeExpressionPtr<ResultType>(fun, args.ptr()...));
        }

        template <class T, class ...ArgTs>
        inline Expression<T> PerformOperation(T fun(ArgTs...), const Expression<ArgTs> & ... args){
            return Expression<T>(MakeExpressionPtr<T>(fun, args.ptr()...));
        }


        // standard operations
        template <class A>
        inline auto operator - (const Expression<A> & a) 
            -> Expression<decltype(-a.eval())> {
            return PerformOperation([](const A& aa){return -aa; }, a);
        }

        template <class A, class B>
        inline auto operator + (const Expression<A> & a, const Expression<B> & b) 
            -> Expression<decltype(a.eval() + b.eval())> {
            return PerformOperation([](const A & aa, const B & bb){return aa + bb; }, a, b);
        }

        template <class A, class B>
        inline auto operator - (const Expression<A> & a, const Expression<B> & b)
            -> Expression<decltype(a.eval() - b.eval())> {
            return PerformOperation([](const A & aa, const B & bb){return aa - bb; }, a, b);
        }

        template <class A, class B>
        inline auto operator * (const Expression<A> & a, const Expression<B> & b)
            -> Expression<decltype(a.eval() * b.eval())> {
            return PerformOperation([](const A & aa, const B & bb){return aa * bb; }, a, b);
        }

        template <class A, class B>
        inline auto operator / (const Expression<A> & a, const Expression<B> & b)
            -> Expression<decltype(a.eval() / b.eval())> {
            return PerformOperation([](const A & aa, const B & bb){return aa / bb; }, a, b);
        }


        inline Expression<long double> sqrt(const Expression<long double> & a) {
            return PerformOperation(std::sqrtl, a);
        }

        inline Expression<float> sqrt(const Expression<float> & a) {
            return PerformOperation(std::sqrtf, a);
        }

        template <class T, int N>
        inline Expression<Vec<T, N>> sqrt(const Expression<Vec<T, N>> & a) {
            return PerformOperation([](const Vec<T, N> & aa){
                Vec<T, N> r;
                for (int i = 0; i < N; i++){
                    r[i] = sqrt(aa[i]);
                }
                return r;
            }, a);
        }

    }
}
 
#endif