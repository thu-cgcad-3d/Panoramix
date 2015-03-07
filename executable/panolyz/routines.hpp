#ifndef PANOLYZ_ROUTINES_HPP
#define PANOLYZ_ROUTINES_HPP

#include "../../src/core/basic_types.hpp"

namespace panolyz {

    enum PanolyzAlgorithm {
        Panorama_v1,
        Normal_v1
    };

    template <PanolyzAlgorithm Algo>
    void Routine(const panoramix::core::Image & image,
        std::integral_constant<PanolyzAlgorithm, Algo> = {});

#define ROUTINE_FOR_ALGORITHM(Algo) \
    template <> \
    void Routine(const panoramix::core::Image & image, \
    std::integral_constant<PanolyzAlgorithm, Algo>)


    // declare routines

    ROUTINE_FOR_ALGORITHM(Panorama_v1);


    ROUTINE_FOR_ALGORITHM(Normal_v1);

}

 
#endif