#pragma once
#ifndef PANORAMIX_DERIV_OPTRAITS_SUM_HPP
#define PANORAMIX_DERIV_OPTRAITS_SUM_HPP

#include "optraits_base.hpp"

namespace panoramix {
    namespace deriv {

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


    }
}
 
#endif