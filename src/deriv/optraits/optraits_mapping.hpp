#pragma once
#ifndef PANORAMIX_DERIV_OPTRAITS_MAPPING_HPP
#define PANORAMIX_DERIV_OPTRAITS_MAPPING_HPP

#include "optraits_base.hpp"

namespace panoramix {
    namespace deriv {

          // mapping
        namespace {
            template <class To, class From, class MappingFunT, class ReverseMappingFunT>
            struct MappingTraits : public OpTraitsBase<To, From> {
                static_assert(IsStorageType<To>::value && IsStorageType<From>::value, "To and From must both be storage types!");
                inline explicit MappingTraits(MappingFunT mf, ReverseMappingFunT rmf, 
                    const std::string & nm = "map", const std::string & rnm = "map") 
                    : mappingFun(mf), reverseMappingFun(rmf), name(nm), reversedName(rnm) {}
                inline To value(ResultType<From> from) const {
                    return mappingFun(from);
                }
                inline void derivatives(
                    Expression<To> to,
                    DerivativeExpression<To> sumOfDOutputs,
                    OriginalAndDerivativeExpression<From> from) const {
                    from.second = ComposeExpression(MappingTraits<From, To, ReverseMappingFunT, MappingFunT>
                        (reverseMappingFun, mappingFun), sumOfDOutputs);
                }
                virtual ostream & toString(ostream & os) const { os << name; return os; }
                MappingFunT mappingFun;
                ReverseMappingFunT reverseMappingFun;
                std::string name, reversedName;
            };
        }

        template <class To, class From, class MappingFunT, class ReverseMappingFunT>
        inline Expression<To> composeMappingFunction(const Expression<From> & from, 
            MappingFunT mf, ReverseMappingFunT rmf, const std::string & nm = "map", const std::string & rnm = "map") {
            return ComposeExpression(MappingTraits<To, From, MappingFunT, ReverseMappingFunT>(mf, rmf, nm, rnm), from);
        }


    }
}
 
#endif