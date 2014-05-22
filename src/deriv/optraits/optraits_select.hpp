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

    }
}
 
#endif