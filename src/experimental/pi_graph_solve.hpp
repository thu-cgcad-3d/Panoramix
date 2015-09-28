#pragma once

#include "../misc/matlab_api.hpp"
#include "pi_graph.hpp"
#include "pi_graph_annotation.hpp"

namespace pano {
    namespace experimental {

        void BuildConstraintGraph(PIGraph & mg);

        void SolvePIGraph(int ccid, PIGraph & mg, misc::Matlab & matlab, int tryNum);

        void ReconstructLayoutAnnotation(PILayoutAnnotation & anno, misc::Matlab & matlab);
      
    }
}