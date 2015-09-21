#pragma once

#include "../misc/matlab_api.hpp"
#include "pi_graph.hpp"

namespace pano {
    namespace experimental {

        void SolvePIGraph(PIGraph & mg, misc::Matlab & matlab);

    }
}