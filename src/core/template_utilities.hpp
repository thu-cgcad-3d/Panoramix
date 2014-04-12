#ifndef PANORAMIX_CORE_TEMPLATE_UTILITIES_HPP
#define PANORAMIX_CORE_TEMPLATE_UTILITIES_HPP
 
#include <tuple>

namespace panoramix {
    namespace core {

        
        // a templated integer sequence
        template<int ...>
        struct Sequence { };

        // use SequenceGenerator<N>::type to deduct Sequence<0, 1, 2, 3, ..., N-1>
        template<int N, int ...S>
        struct SequenceGenerator : SequenceGenerator<N - 1, N - 1, S...> { };
        template<int ...S>
        struct SequenceGenerator<0, S...> {
            typedef Sequence<S...> type;
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