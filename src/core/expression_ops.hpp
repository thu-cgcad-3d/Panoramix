#ifndef PANORAMIX_CORE_EXPRESSION_OPS_HPP
#define PANORAMIX_CORE_EXPRESSION_OPS_HPP

#include "expression.hpp"
 
namespace panoramix {
    namespace core {

        template <class OutputT, class ...InputTs>
        struct OpTraitsBase {
            using OutputType = OutputT;
            static const int InputN = sizeof...(InputTs);
            using InputIdices = typename SequenceGenerator<sizeof...(InputTs)>::type;
            using InputTuple = std::tuple<InputTs...>;
        };

        template <class OpTraitsT, class OutputT, class ...InputTs>
        struct FunctionalOp : public OpWithACache<OutputT> {
        private:
            using InputIdices = typename SequenceGenerator<sizeof...(InputTs)>::type;
            using InputTuple = std::tuple<InputTs...>;

        public:
            inline FunctionalOp(const OutputT & initialValue) : OpWithACache<OutputT>(initialValue) {}

        private:
            template <int ...S>
            inline void evalUsingSequence(ExpressionGraph const * const g, std::vector<EHandle> && inputs, 
                Sequence<S...>) {
                OpTraitsT::evalTrait(cache, g->op(inputs[S]).as<typename std::tuple_element<S, InputTuple>::type>().value()...);
            }

            template <int ...S>
            inline std::vector<EHandle> makeInputDerivativesUsingSequence(ExpressionGraph * const g,
                EHandle self,
                std::vector<EHandle> && inputs,
                EHandle outputDerivsSum, 
                Sequence<S...>) const {
               
                std::tuple<Expression<InputTs>...> inputDExprs; // make a tuple to store resulted input derivatives
                
                OpTraitsT::makeInputDerivativesTrait(
                    Expression<OutputT>(outputDerivsSum, g), // Expression<OutputT>
                    Expression<InputTs>(inputs.at(S), g)..., // Expression<Inputs> ...
                    std::get<S>(inputDExprs)...); // Expression<Inputs> &

                return std::vector<EHandle>{std::get<S>(inputDExprs).handle()...};
            }

        public:
            // get name
            virtual std::ostream & toString(std::ostream & os) const { os << "Func<" << typeid(OpTraitsT).name() << ">"; return os; }

            // evaluate and get result
            virtual void eval(ExpressionGraph const * const g, std::vector<EHandle> && inputs) override {
                assert(inputs.size() == sizeof...(InputTs));
                evalUsingSequence(g, std::move(inputs), InputIdices());
            }

            // make input derivative expressions based on inputs and output derivative expressions
            // returns derivative Ops corresponding to inputs
            virtual std::vector<EHandle> makeInputDerivatives(ExpressionGraph * const g,
                EHandle self,
                std::vector<EHandle> && inputs,
                EHandle outputDerivsSum) const override {

                assert(inputs.size() == sizeof...(InputTs));
                return makeInputDerivativesUsingSequence(g, self, std::move(inputs), outputDerivsSum, InputIdices());
            }

        };

        // compose Expression using FunctionalOp
        template <class OpTraitsT, class InputT, class ...InputTs>
        inline Expression<typename OpTraitsT::OutputType> ComposeExpression(const typename OpTraitsT::OutputType& initialValue, 
            const Expression<InputT> & firstInput,
            const Expression<InputTs> & ... inputs){
            return Expression<typename OpTraitsT::OutputType>(firstInput.g()->
                addNode(std::make_shared<FunctionalOp<OpTraitsT, typename OpTraitsT::OutputType, InputT, InputTs...>>(initialValue), {
                firstInput.handle(), inputs.handle()... 
            }), firstInput.g());
        }




        ///////////////////////////////////////////////////////////////////////////////////////
        //// operators
        ///////////////////////////////////////////////////////////////////////////////////////

        using namespace Eigen;

        template <class T, class K>
        struct ConvertTraits : public OpTraitsBase<T, K> {
            //static_assert(std::is_assignable<T&, const K&>::value, "conversion unavailable!");
            inline static void evalTrait(T & result, const K & from) {
                result = from;
            }
            inline static void makeInputDerivativesTrait(
                Expression<T> dresult,
                Expression<K> from,
                Expression<K> & dfrom) {
                dfrom = ConvertTo<K>(dresult);
            }
        };

        template <class T, class K>
        std::enable_if_t<!std::is_same<T, K>::value, Expression<T>> ConvertTo(const Expression<K> & from) {
            return ComposeExpression<ConvertTraits<T, K>>(from.value(), from);
        }
        template <class T>
        Expression<T> ConvertTo(const Expression<T> & from) { 
            return from; 
        }



        // use the shape of the first input value and fill it with the element value from the second input
        template <class T>
        struct FillWithScalarTraits : public OpTraitsBase<T, T, ElementType<T>> {
            using ScalarType = ElementType<T>;
            inline static void evalTrait(T & result, const T & original, const ScalarType & s) {
                result = original;
                Fill<T>(result, s);
            }
            inline static void makeInputDerivativesTrait(
                Expression<T> dresult,
                Expression<T> original, Expression<ScalarType> s,
                Expression<T> & doriginal, Expression<ScalarType> & ds) {
                doriginal = SetConstant<T>(dresult, 0);
                ds = SumElements<T>(dresult);
            }
        };

        template <class T>
        Expression<T> FillWithScalar(const Expression<T> & f, const Expression<ElementType<T>> & s) {
            T fval = f.value();
            Fill(fval, s.value());
            return ComposeExpression<FillWithScalarTraits<T>>(fval, f, s);
        }



        // sum elements of an input: DataT -> ElementType<DataT>
        template <class T>
        struct SumElementsTraits : public OpTraitsBase<ElementType<T>, T> {
            using ScalarType = ElementType<T>;
            inline static void evalTrait(ScalarType & result, const T & original) {
                GetElementsSum<T>(original, result);
            }
            inline static void makeInputDerivativesTrait(
                Expression<ScalarType> dresult,
                Expression<T> original,
                Expression<T> & doriginal) {
                doriginal = FillWithScalar(original, dresult);
            }
        };

        template <class T>
        Expression<ElementType<T>> SumElements(const Expression<T> & e) {
            ElementType<T> sum;
            GetElementsSum<T>(e.value(), sum);
            return ComposeExpression<SumElementsTraits<T>>(sum, e);
        }


        
        // transpose
        template <class T>
        struct TransposeTraits : public OpTraitsBase<TransposeResultType<T>, T> {
            inline static void evalTrait(TransposeResultType<T> & result, const T & a) {
                Transpose(a, result);
            }
            inline static void makeInputDerivativesTrait(
                Expression<TransposeResultType<T>> dresult,
                Expression<T> a,
                Expression<T> & da) {
                da = Transpose(dresult);
            }
        };

        template <class T>
        inline std::enable_if_t<IsStaticallyScalar<T>::value && std::is_same<T, TransposeResultType<T>>::value, Expression<T>>
            Transpose(const Expression<T> & e) { return e; }

        template <class T>
        inline auto Transpose(const Expression<T> & e) 
            -> Expression<EvaluatedType<decltype(Transpose(std::declval<T>()))>> {
            using ResultType = EvaluatedType<decltype(Transpose(std::declval<T>()))>;
            return ComposeExpression<ResultType>(Transpose(e.value()), e);
        }
        

        // elementwise linear combination trait
        namespace {

            template <class A>
            inline A Sum(const A & a) {
                return a;
            }

            template <class A, class ...Others>
            inline auto Sum(const A & a, const Others &... others) 
                -> decltype(a + Sum(others...)){
                return a + Sum(others...);
            }


        }

        //template <class T, class InputTupleT, int ...Coefs>
        //struct ElementWiseLinearCombinationTraits : public OpTraitsBase<T, InputTs...> {
        //    static_assert(std::tuple_size<InputTuple>::value == sizeof...(Coefs), 
        //    "number of input types and number of coefficients mismatched!");            

        //    inline static void evalTrait(T & result, const InputTs & ... inputs) {
        //        result = Sum(inputs * Coefs ...);
        //    }

        //    inline static void makeInputDerivativesTrait(
        //        Expression<T> dresult,
        //        Expression<InputTs> ... inputs,
        //        Expression<InputTs> & ... dinputs) {
        //        //TODO
        //    }
        //};



        // general elementwise plus2
        template <class T, class T1, class T2>
        struct PlusTraits : public OpTraitsBase<T, T1, T2> {
            inline static void evalTrait(T & result, const T1 & a, const T2 & b) {
                result = a + b;
            }
            inline static void makeInputDerivativesTrait(
                Expression<T> dresult, 
                Expression<T1> a, Expression<T2> b,
                Expression<T1> & da, Expression<T2> & db) {
                da = ConvertTo<T1>(dresult);
                db = ConvertTo<T2>(dresult);
            }
        };

        template <class T1, class T2>
        inline auto operator + (Expression<T1> a, Expression<T2> b) 
            -> Expression<EvaluatedType<decltype(std::declval<T1>() + std::declval<T2>())>> {
            using T = EvaluatedType<decltype(std::declval<T1>() + std::declval<T2>())>;
            return ComposeExpression<PlusTraits<T, T1, T2>>(a.value() + b.value(), a, b);
        }


        // general mult2
        template <class T, class T1, class T2>
        struct MultTraits : public OpTraitsBase<T, T1, T2> {
            inline static void evalTrait(T & result, const T1 & a, const T2 & b) {
                result = a * b; // element wise mult OR matrix mult !!!
            }

            inline static void makeInputDerivativesTrait(
                Expression<T> dresult, 
                Expression<T1> a, Expression<T2> b,
                Expression<T1> & da, Expression<T2> & db) {
                da = ConvertTo<T1>(dresult * Transpose<T2>(b)); // !
                db = ConvertTo<T2>(Transpose<T1>(a) * dresult);
            }
        };

        template <class T1, class T2>
        inline auto operator * (const Expression<T1> & a, const Expression<T2> & b) 
            -> Expression<EvaluatedType<decltype(std::declval<T1>() * std::declval<T2>())>>{
            using T = EvaluatedType<decltype(std::declval<T1>() * std::declval<T2>())>;
            return ComposeExpression<MultTraits<T, T1, T2>>(a.value() * b.value(), a, b);
        }


        // general div
        template <class T, class T1, class T2>
        struct DivTraits : public OpTraitsBase<T, T1, T2> {
            inline static void evalTrait(T & result, const T1 & a, const T2 & b) {
                result = a / b;
            }

            inline static void makeInputDerivativesTrait(
                Expression<T> dresult,
                Expression<T1> a, Expression<T2> b,
                Expression<T1> & da, Expression<T2> & db) {
                da = ConvertTo<T1>(dresult / b);
                //db = ConvertTo<T2>()
            }
        };
        


        // exponent
        template <class T>
        struct ExpTraits : public OpTraitsBase<T, T> {
            inline static void evalTrait(T & output, const T & input) {
                output = exp(input);
            }

            inline static void makeInputDerivativesTrait(
                Expression<T> outputD, 
                Expression<T> input,
                Expression<T> & inputD) {
                inputD = outputD * exp(input);
            }
        };

        template <class T>
        inline Expression<T> exp(const Expression<T> & e) {
            return ComposeExpression<ExpTraits<T>>(exp(e.value()), e);
        }


        // log
        template <class T>
        struct LogTraits : public OpTraitsBase<T, T> {
            inline static void evalTrait(T & output, const T & input) {
                output = log(input);
            }
            inline static void makeInputDerivativesTrait(
                Expression<T> dresult,
                Expression<T> input,
                Expression<T>& dinput) {
                
            }
        };



    }
}
 
#endif