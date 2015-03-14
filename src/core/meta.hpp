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






        // determine whether T is a container
        namespace {
            template <class T>
            struct IsContainerImp {
                template <class TT>
                static auto test(int) -> decltype(
                    std::begin(std::declval<TT>()),
                    std::end(std::declval<TT>()),
                    std::true_type()
                    );
                template <class>
                static std::false_type test(...);
                static const bool value = std::is_same<decltype(test<T>(0)), std::true_type>::value;
            };
        }

        template <class T>
        struct IsContainer : std::integral_constant<bool, IsContainerImp<T>::value> {};

        // determine whether T can accepts ArgTs as arguments
        namespace {
            template <class T, class ... ArgTs>
            struct IsCallableImp {
                template <class TT, class ... AAs>
                static auto test(int) -> decltype(std::declval<TT>()(std::declval<AAs>()...), std::true_type());
                template <class ...>
                static std::false_type test(...);
                static const bool value = std::is_same<decltype(test<T, ArgTs...>(0)), std::true_type>::value;
            };
        }

        template <class T, class ... ArgTs>
        struct IsCallable : std::integral_constant<bool, IsCallableImp<T, ArgTs...>::value> {};


        // function traits
        template <class T>
        struct IsFunction : std::false_type {};

        template <class ResultT, class ... ArgTs>
        struct IsFunction<ResultT(ArgTs...)> : std::true_type{};

        template <class T>
        struct FunctionTraits : public FunctionTraits<decltype(&T::operator())> {};

        template <class ResultT, class ... ArgTs>
        struct FunctionTraits<ResultT(*)(ArgTs...)> {
            using ResultType = ResultT;
            using ArgumentsTupleType = std::tuple<ArgTs...>;
            enum { ArgumentsNum = sizeof...(ArgTs) };
        };

        template <class ClassT, class ResultT, class... ArgTs>
        struct FunctionTraits<ResultT(ClassT::*)(ArgTs...) const> {
            // we specialize for pointers to member function
            using ResultType = ResultT;
            using ArgumentsTupleType = std::tuple<ArgTs...>;
            enum { ArgumentsNum = sizeof...(ArgTs) };
        };


        // iterate over
        namespace {
            template <class FunT, class T>
            inline void IterateOverImp(T && t, FunT && fun, std::false_type){
                fun(std::forward<T>(t));
            }
            template <class FunT, class T>
            inline void IterateOverImp(T && t, FunT && fun, std::true_type){
                for (auto && e : t){
                    IterateOver(e, fun);
                }
            }
        }

        template <class FunT, class T>
        inline void IterateOver(T && t, FunT && fun){
            IterateOverImp(std::forward<T>(t), std::forward<FunT>(fun), std::integral_constant<bool, IsContainer<T>::value>());
        }

        template <class FunT, class T1, class T2>
        inline void IterateOver(std::pair<T1, T2> & p, FunT && fun){
            IterateOver(p.first, fun);
            IterateOver(p.second, fun);
        }

        template <class FunT, class T1, class T2>
        inline void IterateOver(const std::pair<T1, T2> & p, FunT && fun){
            IterateOver(p.first, fun);
            IterateOver(p.second, fun);
        }

        template <class FunT, class T1, class T2>
        inline void IterateOver(std::pair<T1, T2> && p, FunT && fun){
            IterateOver(std::move(p.first), fun);
            IterateOver(std::move(p.second), fun);
        }

        namespace {
            template <class FunT, class ArgT>
            inline bool EvalOneArg(FunT && fun, ArgT && arg) { 
                IterateOver(std::forward<ArgT>(arg), std::forward<FunT>(fun)); return true; 
            }
            template <class FunT, class TupleT, int ... I>
            inline void IterateOverTupleUsingSequence(TupleT && t, FunT && fun, Sequence<I...>){
                bool dummy[] = { EvalOneArg(fun, std::get<I>(t))... };
            }
        }

        template <class FunT, class ...Ts>
        inline void IterateOver(std::tuple<Ts...> & t, FunT && fun){
            IterateOverTupleUsingSequence(t, std::forward<FunT>(fun), SequenceGenerator<sizeof...(Ts)>::type());
        }

        template <class FunT, class ...Ts>
        inline void IterateOver(const std::tuple<Ts...> & t, FunT && fun){
            IterateOverTupleUsingSequence(t, std::forward<FunT>(fun), SequenceGenerator<sizeof...(Ts)>::type());
        }

        template <class FunT, class ...Ts>
        inline void IterateOver(std::tuple<Ts...> && t, FunT && fun){
            IterateOverTupleUsingSequence(std::move(t), std::forward<FunT>(fun), SequenceGenerator<sizeof...(Ts)>::type());
        }







        // invocation with tuple
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



        namespace {
            template <class FunT, class TupleT, int ... Is>
            inline auto TupleMapUsingSequence(FunT && fun, const TupleT & t, Sequence<Is...>) -> decltype(std::make_tuple(fun(std::get<Is>(t)) ...)) {
                return std::make_tuple(fun(std::get<Is>(t)) ...);
            }
            /*template <int I, class FunT, class TupleTupleT>
            inline auto TupleOpEach(FunT && fun, TupleTupleT && ts) -> decltype(fun(std::get<I>(ts)...)) {
                return Invoke(fun, )
            }
            template <class FunT, class ... TupleTs, int ... Is>
            inline auto TupleOpUsingSequence(FunT && fun, const TupleTs & ... ts, Sequence<Is...>){
                return std::make_tuple(TupleOpEach<Is>(fun, ts ...)
            }*/
        }

        template <class FunT, class ... Ts>
        inline auto TupleMap(FunT && fun, const std::tuple<Ts...> & t) -> decltype(TupleMapUsingSequence(std::forward<FunT>(fun), t, SequenceGenerator<sizeof...(Ts)>::type())) {
            return TupleMapUsingSequence(std::forward<FunT>(fun), t, SequenceGenerator<sizeof...(Ts)>::type());
        }




        // some common used functors for std classes
        struct SizeFunctor {
            template <class T>
            inline size_t operator()(const T & t) const{
                return t.size();
            }
        };

        struct ReserveFunctor {
            template <class T>
            inline void operator()(T & t, size_t capacity) const {
                t.reserve(capacity);
            }
        };

        struct PushBackFunctor {
            template <class T, class E>
            inline void operator()(T & t, E && e) const {
                t.push_back(std::forward<E>(e));
            }
        };

        struct ClearFunctor {
            template <class T>
            inline void operator()(T & t) const{
                t.clear();
            }
        };

        template <class StreamerT, class SeperatorT>
        struct GeneralStreamFunctor {
            StreamerT streamer;
            SeperatorT seperator;
            template <class T>
            inline void operator()(T && t) const {
                streamer = (streamer << std::forward<T>(t) << seperator);
            }
        };
        
        template <class StreamerT, class SeperatorT>
        inline GeneralStreamFunctor<StreamerT, std::decay_t<SeperatorT>> MakeStreamFunctor(
            StreamerT && streamer, SeperatorT && sep) {
            return GeneralStreamFunctor<StreamerT, std::decay_t<SeperatorT>>{ 
                std::forward<StreamerT>(streamer), 
                    std::forward<SeperatorT>(sep) 
            };
        }



        // to represent any functor type which returns constant value
        template <class RetT, RetT Val>
        struct StaticConstantFunctor {
            static const RetT value = Val;
            template <class ... ParamTs>
            inline RetT operator()(ParamTs && ... params) const { return value; }
        };

        template <class RetT>
        struct ConstantFunctor {
            inline ConstantFunctor(const RetT & v) : value(v) {}
            inline ConstantFunctor(RetT && v) : value(std::move(v)) {}
            template <class ... ParamTs>
            inline RetT operator()(ParamTs && ... params) const { return value; }
            const RetT value;
        };



    }
}
 
#endif