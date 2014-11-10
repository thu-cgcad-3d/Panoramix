#ifndef PANORAMIX_CORE_META_HPP
#define PANORAMIX_CORE_META_HPP
 
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
        struct IsSequence : no {};

        template <int ... S>
        struct IsSequence<Sequence<S...>> : yes {};

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

        // sequence contains element judger
        template <int E, class SequenceT>
        struct SequenceContainsElement {};

        template <int E>
        struct SequenceContainsElement<E, Sequence<>> {
            static const bool value = false;
        };

        template <int E, int ...S>
        struct SequenceContainsElement<E, Sequence<E, S...>> {
            static const bool value = true;
        };

        template <int E, int N, int ...S>
        struct SequenceContainsElement<E, Sequence<N, S...>> {
            static const bool value = SequenceContainsElement<E, Sequence<S...>>::value;
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


        // is tuple judger
        template <class T>
        struct IsTuple : no {};

        template <class ...T>
        struct IsTuple<std::tuple<T...>> : yes {};


        // tuple contains element type judger
        template <class E, class TupleT>
        struct TupleContainsType {};

        template <class E>
        struct TupleContainsType<E, std::tuple<>> {
            static const bool value = false;
        };

        template <class E, class ...Ts>
        struct TupleContainsType<E, std::tuple<E, Ts...>> {
            static const bool value = true;
        };

        template <class E, class T, class ...Ts>
        struct TupleContainsType<E, std::tuple<T, Ts...>> {
            static const bool value = TupleContainsType<E, std::tuple<Ts...>>::value;
        };

        // find location of tuple element type 
        template <class E, class TupleT>
        struct TypeFirstLocationInTuple {};

        template <class E>
        struct TypeFirstLocationInTuple<E, std::tuple<>> {
            static const int value = -1;
        };

        template <class E, class ...Ts>
        struct TypeFirstLocationInTuple<E, std::tuple<E, Ts...>> {
            static const int value = 0;
        };

        template <class E, class T, class ...Ts>
        struct TypeFirstLocationInTuple<E, std::tuple<T, Ts...>> {
            enum  { _v = TypeFirstLocationInTuple<E, std::tuple<Ts...>>::value };
            static const int value = (_v == -1) ? -1 : (_v + 1);
        };



        namespace {
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