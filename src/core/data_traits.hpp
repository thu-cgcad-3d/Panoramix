#ifndef PANOPTIC_CORE_DATA_TRAITS_HPP
#define PANOPTIC_CORE_DATA_TRAITS_HPP

#include "basic_types.hpp"

namespace panoramix {
    namespace core {
        
        template <class DataT>
        struct DataTraits {
            static DataT one() { return 1; }
            static DataT zero() { return 0; }
        };



    }
}
 
#endif