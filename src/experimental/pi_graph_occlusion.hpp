#include "pi_graph.hpp"

namespace pano {
    namespace experimental {

        void DetectOcclusions(PIGraph & mg);

        void DetectOcclusions2(PIGraph & mg, 
            double minAngleSizeOfLineInTJunction = DegreesToRadians(3),
            double lambdaShrinkForHLineDetectionInTJunction = 0.2,
            double lambdaShrinkForVLineDetectionInTJunction = 0.1,
            double angleSizeForPixelsNearLines = DegreesToRadians(2));

    }
}