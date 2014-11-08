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

        template <int N, int ...S>
        struct Sequence<N, S...> {
            static const int Count = 1 + Sequence<S...>::Count;
            static const int Sum = N + Sequence<S...>::Sum;
            static const int Product = N * Sequence<S...>::Product;
            static const bool All = (Product != 0);
            static const bool Any = (Sum != 0);
        };

        // sequence type judger
        template <class T>
        struct IsSequence : std::false_type {};

        template <int ... S>
        struct IsSequence<Sequence<S...>> : std::true_type {};

        // sequence element getter
        template <int Id, class SequenceT>
        struct SequenceElement {};

        template <int N, int ...S>
        struct SequenceElement<0, Sequence<N, S...>> {
            static const int value = N;
        };

        template <int Id, int N, int ...S>
        struct SequenceElement<Id, Sequence<N, S...>>{
            static const int value = SequenceElement<Id - 1, Sequence<S...>>::value;
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

    }
}
 
#endif