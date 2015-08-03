#pragma once


#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "macros.hpp"

namespace pano {
    namespace core {


        namespace ring_config {

            template <class T>
            struct RadianConfig {
                static T LowerBound() { return 0.0; }
                static T UpperBound() { return M_PI * 2; }
            };

            template <class T>
            struct AngleConfig {
                static T LowerBound() { return 0.0; }
                static T UpperBound() { return 360.0; }
            };

        }


        template <class T, class RepT, class ConfigT>
        class Ring {
            static_assert(std::is_unsigned<RepT>::value, "");
            static_assert(!std::is_same<T, RepT>::value, "");
        public:
            Ring() : _rep(0) {}
            Ring(T v) : _rep(toRep(v)) {}

            operator T() const { return toValue(_rep); }

            Ring & operator += (Ring r) { _rep += r._rep; return *this; }
            Ring & operator -= (Ring r) { _rep -= r._rep; return *this; }
            Ring & operator *= (T v) { _rep *= v; return *this; }
            Ring & operator /= (T v) { _rep /= v; return *this; }

            bool operator == (Ring r) const { return _rep == r._rep; }
            bool operator != (Ring r) const { return _rep != r._rep; }

            Ring fromRep(RepT rep) { Ring r; r._rep = rep; return r; }

            Ring operator + (Ring r) const { return Ring((RepT)(_rep + r._rep)); }
            Ring operator + (T r) const { return Ring((RepT)(_rep + toRep(r))); }
            Ring operator - (Ring r) const { return Ring((RepT)(_rep - r._rep)); }
            Ring operator - (T r) const { return Ring((RepT)(_rep - toRep(r))); }
            Ring operator * (T r) const { return Ring((RepT)(_rep * r)); }
            Ring operator / (T r) const { return Ring((RepT)(_rep / r)); }

            Ring operator -() const { return Ring(-rep); }

            bool operator < (Ring r) const {
                RepT diff = r._rep - _rep;
                return 0 < diff && diff < (std::numeric_limits<RepT>::max() >> 1);
            }
            bool operator <= (Ring r) const {
                RepT diff = r._rep - _rep;
                return diff < (std::numeric_limits<RepT>::max() >> 1);
            }


        public:
            static T toValue(RepT rep) {
                static T P = ConfigT::UpperBound() - ConfigT::LowerBound();
                T ratio = T(rep) / (T(std::numeric_limits<RepT>::max()) + 1);
                return P * ratio +
                    ConfigT::LowerBound();
            }
            static RepT toRep(T value){
                static T P = ConfigT::UpperBound() - ConfigT::LowerBound();
                return RepT((value - ConfigT::LowerBound()) / P *
                    (T(std::numeric_limits<RepT>::max()) + 1));
            }

        private:
            explicit Ring(RepT r) : _rep(r) {}
            RepT _rep;
        };


        using Radian = Ring<double, uint32_t, ring_config::RadianConfig<double>>;
        using Angle = Ring<double, uint32_t, ring_config::AngleConfig<double>>;

    }
}

