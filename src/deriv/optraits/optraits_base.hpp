#pragma once
#ifndef PANORAMIX_DERIV_OPTRAITS_BASE_HPP
#define PANORAMIX_DERIV_OPTRAITS_BASE_HPP

#include "../expression.hpp"
 
namespace panoramix {
    namespace deriv {

        using std::ostream;

        
        // represents std::pair<OriginalExpression, DerivativeExpression &>
        // used in traits by ComposeExpression
        template <class T>
        using OriginalAndDerivativeExpression = std::pair<Expression<T>, DerivativeExpression<T> &>;
        template <class T>
        inline OriginalAndDerivativeExpression<T> MakeOriginalAndDerivative(Expression<T> && input, DerivativeExpression<T> & dinput) {
            return std::pair<Expression<T>, DerivativeExpression<T> &>{input, dinput};
        }

        // represents an array of OriginalExpression and DerivativeExpression pairs
        template <class T>
        using OriginalAndDerivativeExpressionTable = std::vector<std::pair<Expression<T>, DerivativeExpression<T>>>;

        // the base class for all traits used by ComposeExpression
        template <class OutputT, class ...InputTs>
        struct OpTraitsBase {
            using OutputType = OutputT;
            using OutputExpressionType = Expression<OutputT>;

            static const int InputsNumber = sizeof...(InputTs);
            using InputIndices = typename SequenceGenerator<sizeof...(InputTs)>::type;
            using InputTuple = std::tuple<InputTs...>;

            template <int I>
            struct InputTypeStruct {
                using type = typename std::tuple_element<I, InputTuple>::type;
            };
            template <int I>
            using InputType = typename InputTypeStruct<I>::type;

            // to string
            virtual ostream & toString(ostream & os) const { os << "unknown op"; return os; }
        };

        
        namespace {
            // the specific op class used by ComposeExpression
            template <class OpTraitsT>
            class FunctionalOp : public OpBaseType<typename OpTraitsT::OutputType> {
            private:
                using InputIndices = typename OpTraitsT::InputIndices;

                using OutputType = typename OpTraitsT::OutputType;
                using InputTuple = typename OpTraitsT::InputTuple;
                static const int InputsNumber = OpTraitsT::InputsNumber;
                template <int I>
                struct InputTypeStruct {
                    using type = typename std::tuple_element<I, InputTuple>::type;
                };
                template <int I>
                using InputType = typename InputTypeStruct<I>::type;

            public:
                inline explicit FunctionalOp(OpTraitsT t) : OpBaseType<typename OpTraitsT::OutputType>(), _opTraits(t) {}

            private:
                template <int ...S>
                inline OutputType valueUsingSequence(std::vector<EHandle> && inputs, Sequence<S...>) const {
                    assert(inputs.size() == InputsNumber);
                    // inputs -> ResultType<inputs> -> output
                    return _opTraits.value(TraitsAboutOp<InputType<S>>::Result(graph->op(inputs[S])) ...);
                }

                template <int ...S>
                inline std::vector<EHandle> backPropagateGradientUsingSequence(
                    std::vector<EHandle> && inputs,
                    EHandle sumOfDOutputs,
                    Sequence<S...>) const {

                    // inputs + doutputs -> dinputs
                    std::tuple<DerivativeExpression<InputType<S>>...> dinputs;
                    _opTraits.derivatives(
                        graph->as<OutputType>(self), // output
                        graph->asDerived<OutputType>(sumOfDOutputs), // DOutput
                        MakeOriginalAndDerivative(graph->as<InputType<S>>(inputs[S]), std::get<S>(dinputs))...);
                    return std::vector<EHandle>{std::get<S>(dinputs).handle()...};
                }

            public:
                virtual std::ostream & toString(std::ostream & os) const {
                    return _opTraits.toString(os);
                }
                virtual OutputType value() const override {
                    return valueUsingSequence(graph->inputs(self), InputIndices());
                }
                virtual std::vector<EHandle> backPropagateGradient(EHandle sumOfDOutputs) const {
                    return backPropagateGradientUsingSequence(graph->inputs(self), sumOfDOutputs, InputIndices());
                }
                OpTraitsT _opTraits;
            };
        }

        // compose an Expression given traits of the operation
        template <class OpTraitsT, class InputT, class ...InputTs>
        inline Expression<typename OpTraitsT::OutputType> ComposeExpression(
            OpTraitsT opTraits,
            const Expression<InputT> & firstInput,
            const Expression<InputTs> & ... inputs){
            static_assert(OpTraitsT::InputsNumber == 1 + sizeof...(inputs), "inputs number mismatch!");
            return firstInput.g()->as<typename OpTraitsT::OutputType>(firstInput.g()->
                addNode(std::make_shared<FunctionalOp<OpTraitsT>>(opTraits), {
                firstInput.handle(), inputs.handle()...
            }));
        }

        // compose Expression (with no inputs) given traits of the operation
        template <class OpTraitsT>
        inline Expression<typename OpTraitsT::OutputType> ComposeExpressionWithNoInputs(ExpressionGraph & graph,
            OpTraitsT opTraits){
            static_assert(OpTraitsT::InputsNumber == 0, "inputs number mismatch!");
            return graph.as<typename OpTraitsT::OutputType>(graph.addNode(std::make_shared<FunctionalOp<OpTraitsT>>(opTraits)));
        }

        // compose Expression without derivative definition
        template <class FunT, class InputT, class ...InputTs>
        inline auto ComposeExpressionWithoutDerivativeDefinition(
            FunT && fun, 
            const Expression<InputT> & firstInput, 
            const Expression<InputTs> & ... inputs)
            -> Expression<decltype(fun(firstInput.result(), inputs.result()...))> {

            using RetT = decltype(fun(firstInput.result(), inputs.result()...));
            struct TempTraits : public OpTraitsBase<RetT, InputT, InputTs...> {
                inline explicit TempTraits(FunT && ff) : f(std::forward<FunT>(ff)) {}
                inline RetT value(ResultType<InputT> i1, ResultType<InputTs> ... is) const {
                    return f(i1, is...);
                }
                inline void derivatives(
                    Expression<RetT> output,
                    DerivativeExpression<RetT> sumOfDOutputs,
                    OriginalAndDerivativeExpression<InputT> input1,
                    OriginalAndDerivativeExpression<InputTs> ... inputs) const {
                    assert(false && "derivative not defined!");
                }
                FunT f;
            };
            
            return ComposeExpression(TempTraits(std::forward<FunT>(fun)), firstInput, inputs...);
        }


    }
}
 
#endif