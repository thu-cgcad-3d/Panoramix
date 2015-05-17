#ifndef PANORAMIX_CORE_RING_HPP
#define PANORAMIX_CORE_RING_HPP

#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "macros.hpp"

namespace panoramix {
    namespace core {


        namespace ring_config {

            struct RadianConfig {
                static double LowerBound() { return 0.0; }
                static double UpperBound() { return M_PI * 2; }
            };

            struct AngleConfig {
                static double LowerBound() { return 0.0; }
                static double UpperBound() { return 360.0; }
            };

        }


        template <class RepT, class ConfigT>
        class Ring {
            static_assert(std::is_unsigned<RepT>::value, "");

        public:
            Ring() : _rep(0) {}
            Ring(double v) : _rep(toRep(v)) {}

            operator double() const { return toValue(_rep); }

            Ring & operator += (Ring r) { _rep += r._rep; return *this; }
            Ring & operator -= (Ring r) { _rep -= r._rep; return *this; }
            Ring & operator *= (double v) { _rep *= v; return *this; }
            Ring & operator /= (double v) { _rep /= v; return *this; }

            bool operator == (Ring r) const { return _rep == r._rep; }
            bool operator != (Ring r) const { return _rep != r._rep; }

            Ring fromRep(RepT rep) { Ring r; r._rep = rep; return r; }

            Ring operator + (Ring r) const { return Ring((RepT)(_rep + r._rep)); }
            Ring operator + (double r) const { return Ring((RepT)(_rep + toRep(r))); }
            Ring operator - (Ring r) const { return Ring((RepT)(_rep - r._rep)); }
            Ring operator - (double r) const { return Ring((RepT)(_rep - toRep(r))); }
            Ring operator * (double r) const { return Ring((RepT)(_rep * r)); }
            Ring operator / (double r) const { return Ring((RepT)(_rep / r)); }

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
            static double toValue(RepT rep) {
                static double P = ConfigT::UpperBound() - ConfigT::LowerBound();
                double ratio = double(rep) / (double(std::numeric_limits<RepT>::max()) + 1);
                return P * ratio +
                    ConfigT::LowerBound();
            }
            static RepT toRep(double value){
                static double P = ConfigT::UpperBound() - ConfigT::LowerBound();
                return RepT((value - ConfigT::LowerBound()) / P *
                    (double(std::numeric_limits<RepT>::max()) + 1));
            }

        private:
            explicit Ring(RepT r) : _rep(r) {}
            RepT _rep;
        };

        using Radian = Ring<uint32_t, ring_config::RadianConfig>;

    }
}


#endif
