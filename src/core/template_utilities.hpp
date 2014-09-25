#ifndef PANORAMIX_CORE_TEMPLATE_UTILITIES_HPP
#define PANORAMIX_CORE_TEMPLATE_UTILITIES_HPP
 
#include <tuple>

namespace panoramix {
    namespace core {

        using yes = std::true_type;
        using no = std::false_type;

        
        // a templated integer sequence
        template<int ...>
        struct Sequence { };

        template <>
        struct Sequence<> {
            static const int Count = 0;
            static const int Sum = 0;
            static const int Product = 1;
            static const bool All = (Product != 0);
            static const bool Any = (Sum != 0);
        };
        template <int N>
        struct Sequence<N> {
            static const int Count = 1;
            static const int Sum = N;
            static const int Product = N;
            static const bool All = (Product != 0);
            static const bool Any = (Sum != 0);
        };
        template <int N, int ...S>
        struct Sequence<N, S...> {
            static const int Count = 1 + Sequence<S...>::Count;
            static const int Sum = N + Sequence<S...>::Sum;
            static const int Product = N * Sequence<S...>::Product;
            static const bool All = (Product != 0);
            static const bool Any = (Sum != 0);
        };


        // use SequenceGenerator<N>::type to deduct Sequence<0, 1, 2, 3, ..., N-1>
        template<int N, int ...S>
        struct SequenceGenerator : SequenceGenerator<N - 1, N - 1, S...> { };
        template<int ...S>
        struct SequenceGenerator<0, S...> {
            using type = Sequence<S...>;
        };

        // use SequenceRangeGenerator<From, To>::type to deduct Sequence<From, From+1, ..., To-1>
        template <int From, int To, int ...S>
        struct SequenceRangeGenerator : SequenceRangeGenerator<From, To - 1, To - 1, S...>{};
        template <int From, int ...S>
        struct SequenceRangeGenerator<From, From, S...> {
            using type = Sequence<S...>;
        };


        template <class FunctorT, class TupleT, int ...S>
        inline auto InvokeWithEachTupleArg(FunctorT fun, TupleT args, Sequence<S...>)
            -> decltype(fun(std::get<S>(args) ...)) {
            return fun(std::get<S>(args) ...);
        }

        template <class FunctorT, class ElementFunctorT, class TupleT, int ...S>
        inline auto InvokeWithEachTupleArgPreprocessed(FunctorT fun, ElementFunctorT elementPreprocessor, 
            TupleT args, Sequence<S...>)
            -> decltype(fun(std::get<S>(args) ...)) {
            return fun(elementPreprocessor(std::get<S>(args)) ...);
        }


        // invoke a function with tuple args
        template <class FunctorT, class TupleT>
        inline auto Invoke(FunctorT fun, TupleT args)
            -> decltype(InvokeWithEachTupleArg(fun, args,
            typename SequenceGenerator<std::tuple_size<TupleT>::value>::type())) {
            return InvokeWithEachTupleArg(fun, args,
                typename SequenceGenerator<std::tuple_size<TupleT>::value>::type());
        }

        template <class FunctorT, class ElementFunctorT, class TupleT>
        inline auto InvokeWithPreprocessor(FunctorT fun, ElementFunctorT elementPreprocessor, TupleT args)
            -> decltype(InvokeWithEachTupleArgPreprocessed(fun, elementPreprocessor, args,
            typename SequenceGenerator<std::tuple_size<TupleT>::value>::type())) {
            return InvokeWithEachTupleArg(fun, args,
                typename SequenceGenerator<std::tuple_size<TupleT>::value>::type());
        }


        // invoke a function with tuple args without the original return value
        template <class FunctorT, class TupleT>
        inline void InvokeWithoutReturn(FunctorT fun, TupleT args) {
            InvokeWithEachTupleArg(fun, args,
                typename SequenceGenerator<std::tuple_size<TupleT>::value>::type());
        }


    }
}
 
#endif