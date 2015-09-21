#pragma once

#include "pi_graph.hpp"

namespace pano {
    namespace experimental {

        void AttachPrincipleDirectionConstraints(PIGraph & mg);
        void AttachWallConstraints(PIGraph & mg, double angle);

        void DisableTopSeg(PIGraph & mg);
        void DisableBottomSeg(PIGraph & mg);
        void DisableInvalidConstraints(PIGraph & mg);

        void AttachGCConstraints(PIGraph & mg, const Image5d & gc);

 
    }
}