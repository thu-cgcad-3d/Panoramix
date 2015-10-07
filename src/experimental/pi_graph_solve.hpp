#pragma once

#include "../misc/matlab_api.hpp"
#include "pi_graph.hpp"
#include "pi_graph_annotation.hpp"

namespace pano {
    namespace experimental {

#if 0
        void BuildConstraintGraph(PIGraph & mg);


        std::vector<double> InverseDepthCoefficientsOfSegAtDirection(const PIGraph & mg,
            int seg, const Vec3 & direction);
        std::vector<double> InverseDepthCoefficientsOfLineAtDirection(const PIGraph & mg,
            int line, const Vec3 & direction);


        void SolvePIGraph(int ccid, PIGraph & mg, misc::Matlab & matlab, int tryNum);

        void ReconstructLayoutAnnotation(PILayoutAnnotation & anno, misc::Matlab & matlab);
#endif
    }
}