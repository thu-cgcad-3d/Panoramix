#pragma once

#include "pi_graph.hpp"

namespace pano {
    namespace experimental {

        void AttachPrincipleDirectionConstraints(PIGraph & mg, double angle);
        void AttachWallConstraints(PIGraph & mg, double angle, int vertVPId);

        void DisableTopSeg(PIGraph & mg);
        void DisableBottomSeg(PIGraph & mg);
        void DisableInvalidConstraints(PIGraph & mg);

        void AttachGCConstraints(PIGraph & mg, const Image5d & gc);

        void DetectAndApplyOcclusions(PIGraph & mg);

    }
}