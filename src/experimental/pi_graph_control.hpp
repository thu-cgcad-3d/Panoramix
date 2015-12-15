#pragma once

#include "pi_graph.hpp"

namespace pano {
namespace experimental {

void AttachPrincipleDirectionConstraints(PIGraph<PanoramicCamera> &mg);
void AttachWallConstraints(PIGraph<PanoramicCamera> &mg, double angle);

void DisableTopSeg(PIGraph<PanoramicCamera> &mg);
void DisableBottomSeg(PIGraph<PanoramicCamera> &mg);

void AttachGCConstraints(PIGraph<PanoramicCamera> &mg, const Image5d &gc,
                         double clutterThres = 0.7, double wallThres = 0.5,
                         bool onlyConsiderBottomHalf = true);
}
}