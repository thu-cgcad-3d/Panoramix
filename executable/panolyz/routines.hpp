#ifndef PANOLYZ_ROUTINES_HPP
#define PANOLYZ_ROUTINES_HPP

#include "../../src/core/basic_types.hpp"

namespace panolyz {

    enum PanolyzAlgorithm {
        PanoramaIndoor_v1,
        NormalIndoor_v1,
        PanoramaOutdoor_v1,
        NormalOutdoor_v1
    };

    template <PanolyzAlgorithm Algo>
    void Routine(const panoramix::core::Image & image, const std::string & tag, 
        std::integral_constant<PanolyzAlgorithm, Algo> = {});

#define ROUTINE_FOR_ALGORITHM(Algo) \
    template <> \
    void Routine(const panoramix::core::Image & image, const std::string & tag, \
    std::integral_constant<PanolyzAlgorithm, Algo>)


    // declare routines

    ROUTINE_FOR_ALGORITHM(PanoramaIndoor_v1);
    /*ROUTINE_FOR_ALGORITHM(NormalIndoor_v1);
    ROUTINE_FOR_ALGORITHM(PanoramaOutdoor_v1);
    ROUTINE_FOR_ALGORITHM(NormalOutdoor_v1);*/

}

 
#endif