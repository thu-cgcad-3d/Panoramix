#ifndef PANORAMIX_CORE_RING_HPP
#define PANORAMIX_CORE_RING_HPP

#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace panoramix {
    namespace core {

        namespace ring_config {

            struct RadianConfig {
                __forceinline static double LowerBound() { return 0.0; }
                __forceinline static double UpperBound() { return M_PI * 2; }
            };

            struct AngleConfig {
                __forceinline static double LowerBound() { return 0.0; }
                __forceinline static double UpperBound() { return 360.0; }
            };

        }


        template <class RepT, class ConfigT>
        class Ring {
            static_assert(std::is_unsigned<RepT>::value, "");

        public:
            __forceinline Ring() : _rep(0) {}
            __forceinline Ring(double v) : _rep(toRep(v)) {}
            
            __forceinline operator double() const { return toValue(_rep); }

            __forceinline Ring & operator += (Ring r) { _rep += r._rep; return *this; }
            __forceinline Ring & operator -= (Ring r) { _rep -= r._rep; return *this; }
            __forceinline Ring & operator *= (double v) { _rep *= v; return *this; }
            __forceinline Ring & operator /= (double v) { _rep /= v; return *this; }

            __forceinline bool operator == (Ring r) const { return _rep == r._rep; }
            __forceinline bool operator != (Ring r) const { return _rep != r._rep; }

            __forceinline Ring fromRep(RepT rep) { Ring r; r._rep = rep; return r; }

            __forceinline Ring operator + (Ring r) const { return Ring((RepT)(_rep + r._rep)); }
            __forceinline Ring operator + (double r) const { return Ring((RepT)(_rep + toRep(r))); }
            __forceinline Ring operator - (Ring r) const { return Ring((RepT)(_rep - r._rep)); }
            __forceinline Ring operator - (double r) const { return Ring((RepT)(_rep - toRep(r))); }
            __forceinline Ring operator * (double r) const { return Ring((RepT)(_rep * r)); }
            __forceinline Ring operator / (double r) const { return Ring((RepT)(_rep / r)); }

            __forceinline Ring operator -() const { return Ring(-rep); }


        private:
            __forceinline static double toValue(RepT rep) {
                static double P = ConfigT::UpperBound() - ConfigT::LowerBound();
                double ratio = double(rep) / (double(std::numeric_limits<RepT>::max()) + 1);
                return P * ratio  +
                    ConfigT::LowerBound();
            }
            __forceinline static RepT toRep(double value){
                static double P = ConfigT::UpperBound() - ConfigT::LowerBound();
                return RepT((value - ConfigT::LowerBound()) / P *
                    (double(std::numeric_limits<RepT>::max()) + 1));
            }

        private:
            __forceinline explicit Ring(RepT r) : _rep(r) {}
            RepT _rep;
        };

        using Radian = Ring<uint32_t, ring_config::RadianConfig>;

    }
}


#endif
