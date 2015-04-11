#ifndef PANORAMIX_MISC_TOOLS_HPP
#define PANORAMIX_MISC_TOOLS_HPP

#include <algorithm>
#include "../core/basic_types.hpp"

namespace panoramix {
    namespace misc {



        
        namespace dataset {

            namespace YorkUrbanDB {
               

                struct GroundTruth {
                    //std::vector<core::Line2> lines;
                    //std::vector<int> vp_association;
                    std::vector<core::Vec3> vp_orthogonal;
                };

                GroundTruth LoadGroundTruth(const std::string & camParams);
            }



        }
    }
}

#endif