#pragma once
#ifndef PANORAMIX_DERIV_OPTRAITS_COMPOSE_HPP
#define PANORAMIX_DERIV_OPTRAITS_COMPOSE_HPP

#include "optraits_base.hpp"
 
namespace panoramix {
    namespace deriv {

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
        inline auto composeFunction(const std::string & name, ExpressionGraph & graph, FunctorT && fun) -> Expression<decltype(fun())> {
            return ComposeExpressionWithNoInputs(graph, 
                FunctorComposerTraits<FunctorT, decltype(fun())>(name, std::forward<FunctorT>(fun)));
        }


    }
}
 
#endif