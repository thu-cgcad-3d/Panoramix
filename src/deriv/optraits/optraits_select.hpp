#pragma once
#ifndef PANORAMIX_DERIV_OPTRAITS_SELECT_HPP
#define PANORAMIX_DERIV_OPTRAITS_SELECT_HPP
 
#include "optraits_base.hpp"

namespace panoramix {
    namespace deriv {

         // cwiseSelect
        namespace {

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectTraits 
                : public OpTraitsBase<decltype(common::CWiseSelect(
                std::declval<ResultType<CondT>>(),
                std::declval<ResultType<IfT>>(),
                std::declval<ResultType<ElseT>>())), CondT, IfT, ElseT> {

                inline OutputType value(ResultType<CondT> cond, ResultType<IfT> ifval, ResultType<ElseT> elseval) const {
                    return common::CWiseSelect(cond, ifval, elseval);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<CondT> cond,
                    OriginalAndDerivativeExpression<IfT> ifval,
                    OriginalAndDerivativeExpression<ElseT> elseval) const {
                    ifval.second = cwiseSelect(cond.first, sumOfDOutputs, 0.0).eval();
                    elseval.second = cwiseSelect(cond.first, 0.0, sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "cwiseSelect"; return os; }
            };

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectWhenIfRetIsConstTraits 
                : public OpTraitsBase<decltype(common::CWiseSelect(
                std::declval<ResultType<CondT>>(),
                std::declval<ResultType<DataStorageType<IfT>>>(),
                std::declval<ResultType<ElseT>>())), CondT, ElseT> {

                inline explicit CWiseSelectWhenIfRetIsConstTraits(const IfT & ifv) : ifval(ifv) {}
                inline OutputType value(ResultType<CondT> cond, ResultType<ElseT> elseval) const {
                    return common::CWiseSelect(cond, ifval, elseval);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<CondT> cond,
                    OriginalAndDerivativeExpression<ElseT> elseval) const {
                    elseval.second = cwiseSelect(cond.first, 0.0, sumOfDOutputs).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "cwiseSelect[ifthen:" << ifval << "]"; return os; }
                DataStorageType<IfT> ifval;
            };

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectWhenElseRetIsConstTraits 
                : public OpTraitsBase<decltype(common::CWiseSelect(
                std::declval<ResultType<CondT>>(),
                std::declval<ResultType<IfT>>(),
                std::declval<ResultType<DataStorageType<ElseT>>>())), CondT, IfT> {

                inline explicit CWiseSelectWhenElseRetIsConstTraits(const ElseT & elsev) : elseval(elsev) {}
                inline OutputType value(ResultType<CondT> cond, ResultType<IfT> ifval) const {
                    return common::CWiseSelect(cond, ifval, elseval);
                }
                inline void derivatives(
                    Expression<OutputType> output,
                    DerivativeExpression<OutputType> sumOfDOutputs,
                    OriginalAndDerivativeExpression<CondT> cond,
                    OriginalAndDerivativeExpression<IfT> ifval) const {
                    ifval.second = cwiseSelect(cond.first, sumOfDOutputs, 0.0).eval();
                }
                virtual ostream & toString(ostream & os) const { os << "cwiseSelect[elsethen:" << elseval << "]"; return os; }
                DataStorageType<ElseT> elseval;
            };

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectComposer {
                inline static Expression<IfT> Compose(const CondT & cond, const IfT & ifval, const ElseT & elseval) {
                    SHOULD_NEVER_BE_INSTANCIATED();
                }
            };

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectComposer<Expression<CondT>, IfT, Expression<ElseT>> {
                inline static auto Compose(const Expression<CondT> & cond, const IfT & ifval, const Expression<ElseT> & elseval)
                -> decltype(ComposeExpression(CWiseSelectWhenIfRetIsConstTraits<CondT, IfT, ElseT>(ifval), cond, elseval)) {
                    return ComposeExpression(CWiseSelectWhenIfRetIsConstTraits<CondT, IfT, ElseT>(ifval), cond, elseval);
                }
            };

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectComposer<Expression<CondT>, Expression<IfT>, ElseT> {
                inline static auto Compose(const Expression<CondT> & cond, const Expression<IfT> & ifval, const ElseT & elseval)
                -> decltype(ComposeExpression(CWiseSelectWhenElseRetIsConstTraits<CondT, IfT, ElseT>(elseval), cond, ifval)) {
                    return ComposeExpression(CWiseSelectWhenElseRetIsConstTraits<CondT, IfT, ElseT>(elseval), cond, ifval);
                }
            };

            template <class CondT, class IfT, class ElseT>
            struct CWiseSelectComposer<Expression<CondT>, Expression<IfT>, Expression<ElseT>> {
                inline static auto Compose(const Expression<CondT> & cond, const Expression<IfT> & ifval, const Expression<ElseT> & elseval)
                -> decltype(ComposeExpression(CWiseSelectTraits<CondT, IfT, ElseT>(), cond, ifval, elseval)) {
                    return ComposeExpression(CWiseSelectTraits<CondT, IfT, ElseT>(), cond, ifval, elseval);
                }
            };

        }

        template <class CondT, class IfT, class ElseT>
        inline auto cwiseSelect(const CondT & cond, const IfT & ifval, const ElseT & elseval)
            -> decltype(CWiseSelectComposer<CondT, IfT, ElseT>::Compose(cond, ifval, elseval)) {
            return CWiseSelectComposer<CondT, IfT, ElseT>::Compose(cond, ifval, elseval);
        }



        template <class ExprIteratorT>
        inline Expression<typename std::iterator_traits<ExprIteratorT>::value_type::Type> 
            minInRange(ExprIteratorT begin, ExprIteratorT end){
            using T = typename std::iterator_traits<ExprIteratorT>::value_type::Type;
            static_assert(IsStorageType<T>::value, "T must be a storage type");
            struct PassIfMinElementIsAssigned : public deriv::OpBaseType < T > {
                inline PassIfMinElementIsAssigned(const std::vector<EHandle> & in, int aid) : inputs(in), assignedId(aid) {}
                virtual std::ostream & toString(std::ostream & os) const {
                    os << "forwardIfMatch";
                    return os;
                }
                virtual T value() const override {
                    auto curMinIt = std::min_element(inputs.begin(), inputs.end(), [this](EHandle a, EHandle b)->bool {
                        return graph->as<T>(a).result() < graph->as<T>(b).result();
                    });
                    int pos = std::distance(inputs.begin(), curMinIt);
                    std::vector<EHandle> inp = graph->inputs(self);
                    assert(inp.size() == 1);
                    if (pos == assignedId){
                        return graph->as<T>(inp.front()).result();
                    }
                    T v;
                    return common::FillWithScalar(v, 0.0);
                }
                virtual std::vector<EHandle> backPropagateGradient(EHandle sumOfDOutputs) const {
                    assert(false && "not implemented yet");
                    return std::vector<EHandle>();
                }
                std::vector<EHandle> inputs;
                int assignedId;
            };
            struct MinElementOp : public deriv::OpBaseType<T> {
                virtual std::ostream & toString(std::ostream & os) const {
                    os << "minElement";
                    return os;
                }
                virtual T value() const override {
                    std::vector<EHandle> inp = graph->inputs(self);
                    assert(!inp.empty() && "wrong inputs number");
                    //return TraitsAboutOp<From>::Result(graph->op(inp.front()));
                    return graph->as<T>(*std::min_element(inp.begin(), inp.end(), [this](EHandle a, EHandle b)->bool {
                        return graph->as<T>(a).result() < graph->as<T>(b).result();
                    })).result();
                }
                virtual std::vector<EHandle> backPropagateGradient(EHandle sumOfDOutputs) const {
                    std::vector<EHandle> inp = graph->inputs(self);
                    assert(!inp.empty() && "wrong inputs number");
                    std::vector<EHandle> derivs(inp.size());
                    for (auto it = derivs.begin(); it != derivs.end(); ++it){
                        *it = graph->addNode(
                            std::make_shared<PassIfMinElementIsAssigned>(inp, std::distance(derivs.begin(), it)),
                            { sumOfDOutputs });
                    }
                    return derivs;
                }
            };

            std::vector<EHandle> inputHandles;
            for (auto exprIt = begin; exprIt != end; ++exprIt){
                inputHandles.push_back(exprIt->handle());
            }
            return (*begin).g()->as<T>((*begin).g()->addNode(std::make_shared<MinElementOp>(), inputHandles));
        }

    }
}
 
#endif